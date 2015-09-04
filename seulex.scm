;;
;;
;; Chicken Scheme interface to the SEULEX numerical solver for systems
;; of first order ordinary differential equations MY'=F(X,Y).
;; 
;;  The SEULEX code is part of the book
;;
;;        E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
;;         EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
;;        SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS 14,
;;        SPRINGER-VERLAG 1991, SECOND EDITION 1996.
;;
;;  Chicken Scheme code Copyright 2011-2012 Ivan Raikov.
;;
;; 
;;  Redistribution and use in source and binary forms, with or without
;;  modification, are permitted provided that the following conditions
;;  are met:
;; 
;;  - Redistributions of source code must retain the above copyright
;;  notice, this list of conditions and the following disclaimer.
;; 
;;  - Redistributions in binary form must reproduce the above
;;  copyright notice, this list of conditions and the following
;;  disclaimer in the documentation and/or other materials provided
;;  with the distribution.
;; 
;;  - Neither name of the copyright holders nor the names of its
;;  contributors may be used to endorse or promote products derived
;;  from this software without specific prior written permission.
;; 
;;  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND THE
;;  CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
;;  INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
;;  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
;;  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR THE
;;  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
;;  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
;;  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
;;  USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
;;  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
;;  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
;;  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
;;  POSSIBILITY OF SUCH DAMAGE.
;;
;;  
;;


(module seulex

	( 
	 
	 seulex-create-solver ;; seulex-create-solver/unsafe
         seulex-destroy-solver seulex-solve seulex-yy
	 seulex-nfcn seulex-njac seulex-nstep seulex-naccpt seulex-nrejct 
	 seulex-ndec seulex-nsol

	 pointer-f64-ref pointer-f64-set! pointer+-f64 
	 seulex_fcn_cb seulex_jac_cb
	 )

	(import scheme chicken foreign)
	(require-extension srfi-4)
	(require-library lolevel srfi-1 )
	(import (only lolevel move-memory! object-evict object-release pointer+ pointer-f64-ref pointer-f64-set! pointer-s32-ref)
		(only srfi-1 fold)
		)


(foreign-declare #<<EOF

/*
     SUBROUTINE SEULEX(N,FCN,IFCN,X,Y,XEND,H,
                       RTOL,ATOL,ITOL,
                       JAC ,IJAC,MLJAC,MUJAC,
                       MAS ,IMAS,MLMAS,MUMAS,
                       SOLOUT,IOUT,
                       WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
*/

typedef void (*rhs_fn) (int*,double*,double*,double*,double*,int*);
typedef void (*jac_fn) (int*,double*,double*,double*,int*,double*,int*);
typedef void (*mas_fn) (int*,double*,int*,double*,int*);
typedef void (*solout_fn) (int*,int*,double*,double*,double*,int*,int*,int*,int*,int*,double*,int*,int*);

extern void seulex_ (int*,rhs_fn,int*,double*,double*,double*,double*,
                     double*,double*,int*,
                     jac_fn,int*,int*,int*,
		     mas_fn,int*,int*,int*,
                     solout_fn,int*,
                     double*,int*,int*,int*,double*,int*,int*);
		     

EOF
)

#>

#include <stdlib.h>
#include <math.h>
#include <assert.h>


static void chicken_panic (C_char *) C_noret;
static void chicken_panic (C_char *msg)
{
  C_word *a = C_alloc (C_SIZEOF_STRING (strlen (msg)));
  C_word scmmsg = C_string2 (&a, msg);
  C_halt (scmmsg);
  exit (5); /* should never get here */
}

static void chicken_throw_exception(C_word value) C_noret;
static void chicken_throw_exception(C_word value)
{
  char *aborthook = C_text("\003sysabort");

  C_word *a = C_alloc(C_SIZEOF_STRING(strlen(aborthook)));
  C_word abort = C_intern2(&a, aborthook);

  abort = C_block_item(abort, 0);
  if (C_immediatep(abort))
    chicken_panic(C_text("`##sys#abort' is not defined"));

  C_save(value);
  C_do_apply(1, abort, C_SCHEME_UNDEFINED);
}

void chicken_error (char *msg, C_word obj) 
{
  size_t msglen;
  C_word *a;
  C_word scmmsg;
  C_word list;

  msglen = strlen (msg);
  a = C_alloc (C_SIZEOF_STRING (msglen) + C_SIZEOF_LIST(2));
  scmmsg = C_string2 (&a, (char *) msg);
  list = C_list(&a, 2, scmmsg, obj);
  chicken_throw_exception(list);
}

typedef struct SEULEX_Solver_Handle_struct
{
    double x;
    int n;

    C_word syy;

    double *yy;
    double abstol;
    double reltol;

    double *work;
    int lwork;

    int *iwork;
    int liwork;

    int ijac;
    int ifcn;
    unsigned int data_index;
    
} SEULEX_Solver_Handle;


SEULEX_Solver_Handle* seulex_create_solver(

  double x_start,

  int variable_number,
  C_word variables, 

  int density_component_number,
  C_word density_components,

  double abstol,
  double reltol,

  int ijac,
  int ifcn,

  unsigned int data_index
)
{
    int flag = 0, i;
    int N, KM, KM2, NRDENS;

    SEULEX_Solver_Handle* solver_handle;
    assert ((solver_handle = malloc (sizeof(struct SEULEX_Solver_Handle_struct))) != NULL);

    solver_handle->x = x_start;

    N = variable_number;
    solver_handle->n = N;

    solver_handle->syy = variables;

    if ((solver_handle->yy = malloc ((1+N)*sizeof(double))) == NULL)
    {
       chicken_error("could not allocate memory for vector yy", C_SCHEME_UNDEFINED);
    }
    memset ((solver_handle->yy),0,(1+N)*sizeof(double));
    memcpy (solver_handle->yy, C_c_f64vector(variables), N*sizeof(double));

    solver_handle->abstol = abstol;
    solver_handle->reltol = reltol;

    NRDENS = density_component_number;
    KM     = 12;
    KM2    = 2+KM*(KM+3)/2;

    solver_handle->lwork  =  N*(2*N+KM+8)+4*KM+20+KM2*NRDENS;
    solver_handle->liwork =  2*N+KM+20+NRDENS;

    if ((solver_handle->work = malloc ((1+solver_handle->lwork)*sizeof(double))) == NULL)
    {
       chicken_error("could not allocate memory for vector work", C_SCHEME_UNDEFINED);
    }
    memset (solver_handle->work, 0, (1+solver_handle->lwork*sizeof(double)));

    if ((solver_handle->iwork = malloc ((1+solver_handle->liwork)*sizeof(int))) == NULL)
    {
       chicken_error("could not allocate memory for vector iwork", C_SCHEME_UNDEFINED);
    }
    memset (solver_handle->iwork, 0, (1+solver_handle->liwork)*sizeof(int));
    
    // number of components for which dense output is required
    solver_handle->iwork[7] = NRDENS;
    if (NRDENS > 0)
    {
       int *dens = C_c_s32vector(density_components);
       for (i=0; i < NRDENS;  i++)
       {
          solver_handle->iwork[i+21] = dens[i];
       }
    }

    solver_handle->ifcn = ifcn;
    solver_handle->ijac = ijac;
    solver_handle->data_index = data_index;

    return solver_handle;
}


void seulex_destroy_solver (SEULEX_Solver_Handle* solver_handle)
{
    free(solver_handle->yy);
    free(solver_handle->work);
    free(solver_handle->iwork);
    free(solver_handle);
    return;
}


<#

(define-foreign-type SEULEXSolverHandle "SEULEX_Solver_Handle")

(define seulex-data-global (make-parameter '( #f )))

(define seulex-fcn-global (make-parameter (lambda _ (begin))))

(define seulex-jac-global (make-parameter (lambda _ (begin))))

;; (N,X,Y,F,RPAR,IPAR)

(define-external (seulex_fcn_cb ((c-pointer int)    n)
				((c-pointer double) x) 
				((c-pointer double) yy) 
				((c-pointer double) rr) 
				((c-pointer double) rpar)
				((c-pointer int)    ipar)) void

  (let ((data-index (pointer-s32-ref ipar)))

    (if (zero? data-index)

	((seulex-fcn-global)  
	 (pointer-f64-ref x) yy rr)

	((seulex-fcn-global) 
	 (pointer-f64-ref x)
	 yy rr (vector-ref (seulex-data-global) data-index ))

	)))

;; (N,X,Y,DFY,LDFY,RPAR,IPAR)

(define-external (seulex_jac_cb ((c-pointer int)    n)
				((c-pointer double) x) 
				((c-pointer double) yy) 
				((c-pointer double) dfy) 
				((c-pointer int)    ldfy) 
				((c-pointer double) rpar)
				((c-pointer int)    ipar)) void

  (let ((data-index (pointer-s32-ref ipar)))

    (if (zero? data-index)

	((seulex-jac-global)  
	 (pointer-f64-ref x) yy dfy)

	((seulex-jac-global) 
	 (pointer-f64-ref x)
	 yy dfy (vector-ref (seulex-data-global) data-index ))

	)))

(define f:seulex
  (foreign-safe-lambda* int (((nonnull-c-pointer SEULEXSolverHandle) solver_handle) 
			     (double tout)) #<<EOF

    int flag;
    double h;
    jac_fn jac;
    int itol, ijac, mujac, imas, mumas, iout;

    flag = 0; h = 1e-3; itol = 0; imas = 0; iout = 0;

    if (solver_handle->ijac > 0)
    {
      jac = seulex_jac_cb;
    } else 
    {
      jac = NULL;
    }

    seulex_ (&(solver_handle->n),
	     seulex_fcn_cb,
	     &(solver_handle->ifcn),
	     &(solver_handle->x),
	     solver_handle->yy,
	     &tout,
	     &h,
	     &(solver_handle->reltol),
	     &(solver_handle->abstol),
	     &itol,
	     jac, &(solver_handle->ijac), &(solver_handle->n), &mujac,
	     NULL, &imas, &(solver_handle->n), &mumas,
	     NULL, &iout, 
	     solver_handle->work, &(solver_handle->lwork), 
	     solver_handle->iwork, &(solver_handle->liwork), 
	     NULL, &(solver_handle->data_index),
	     &flag);

    C_return(flag);
EOF
))


  
(define c-seulex-create-solver
  (foreign-safe-lambda (nonnull-c-pointer SEULEXSolverHandle) 
		  "seulex_create_solver" 
		  double    ;; x_start
		  
		  int ;; variable_number
		  scheme-object ;; variables

		  int ;; density_component_number
		  scheme-object ;; density_components

		  double ;; abstol
		  double ;; reltol

		  int ;; ijac
		  int ;; ifcn
		  unsigned-int ;; user data

		  ))


(define (seulex-create-solver xinit yinit fcn 
	  #!key 
	  (jacobian #f)
	  (autonomous #t)
	  (user-data #f)
	  (density-components (make-s32vector 0))
	  (reltol  1.0e-6)
	  (abstol  1.0e-6)
	  )
  
  (let ((n (f64vector-length yinit))
	(data-index (if user-data (+ 1 (fold max 0 (map car (seulex-data-global)))) 0)))

    (if user-data (seulex-data-global 
		   (list->vector
		    (append (vector->list (seulex-data-global)) 
			    (list user-data)))))

    (let ((fcn1
	   (if user-data

	       (lambda (x yy rr data)  
		 (let ((yy1 (make-f64vector n)))
		   (move-memory! yy yy1 (fx* 8 n))
		   (let ((v (fcn x yy1 data)))
		     (move-memory! v rr (fx* 8 n))
		     )))

	       (lambda (x yy rr)  
		 (let ((yy1 (make-f64vector n 0.0)))
		   (move-memory! yy yy1 (fx* 8 n))
		   (let ((v (fcn x yy1)))
		     (move-memory! v rr (fx* 8 n))
		     )))
	       ))
	  (jac
	   (and jacobian
		(if user-data
		    
		    (lambda (x yy dfy data)  
		      (let ((yy1 (make-f64vector n)))
			(move-memory! yy yy1 (fx* 8 n))
			(let ((v (jacobian x yy1 data)))
			  (move-memory! v dfy (fx* 8 (fx* n n)))
			  )))
		    
		    (lambda (x yy dfy)
		      (let ((yy1 (make-f64vector n 0.0)))
			(move-memory! yy yy1 (fx* 8 n))
			(let ((v (jacobian x yy1)))
			  (move-memory! v dfy (fx* 8 (fx* n n)))
			  )))
		    )))
		  
	  )

      (seulex-fcn-global fcn1)
      (if jac (seulex-jac-global jac))

      (c-seulex-create-solver
       xinit  
       n (object-evict yinit) 
       (s32vector-length density-components) 
       (object-evict density-components)
       abstol reltol
       (if jac 1 0)
       (if autonomous 0 1) 
       data-index)

      )))
		       
			   

(define c-seulex-destroy-solver
  (foreign-safe-lambda void "seulex_destroy_solver" (nonnull-c-pointer SEULEXSolverHandle) ))


(define (seulex-destroy-solver solver)
  (object-release (seulex-syy solver))
  (c-seulex-destroy-solver solver))


(define (seulex-solve handle tout)
  (f:seulex handle tout))


(define seulex-nfcn
  (foreign-safe-lambda* int (((nonnull-c-pointer SEULEXSolverHandle) handle))
#<<EOF
    int result;
    result = handle->iwork[14];
    C_return (result);
EOF
))


(define seulex-njac
  (foreign-safe-lambda* int (((nonnull-c-pointer SEULEXSolverHandle) handle))
#<<EOF
    int result;
    result = handle->iwork[15];
    C_return (result);
EOF
))


(define seulex-nstep
  (foreign-safe-lambda* int (((nonnull-c-pointer SEULEXSolverHandle) handle))
#<<EOF
    int result;
    result = handle->iwork[16];
    C_return (result);
EOF
))


(define seulex-naccpt
  (foreign-safe-lambda* int (((nonnull-c-pointer SEULEXSolverHandle) handle))
#<<EOF
    int result;
    result = handle->iwork[17];
    C_return (result);
EOF
))


(define seulex-nrejct
  (foreign-safe-lambda* int (((nonnull-c-pointer SEULEXSolverHandle) handle))
#<<EOF
    int result;
    result = handle->iwork[18];
    C_return (result);
EOF
))


(define seulex-ndec
  (foreign-safe-lambda* int (((nonnull-c-pointer SEULEXSolverHandle) handle))
#<<EOF
    int result;
    result = handle->iwork[19];
    C_return (result);
EOF
))


(define seulex-nsol
  (foreign-safe-lambda* int (((nonnull-c-pointer SEULEXSolverHandle) handle))
#<<EOF
    int result;
    result = handle->iwork[20];
    C_return (result);
EOF
))

(define c_seulex_yy 
  (foreign-safe-lambda* void (((nonnull-c-pointer SEULEXSolverHandle) handle) 
			      (f64vector result) 
			      (int n))
#<<EOF
   memcpy(result,
	  handle->yy,
	  n*sizeof(double));
EOF
))

(define c_seulex_yy_length
  (foreign-safe-lambda* int (((nonnull-c-pointer SEULEXSolverHandle) handle))
#<<EOF
   C_return (handle->n);
EOF
))

(define (seulex-yy handle)
  (let* ((n (c_seulex_yy_length handle))
	 (v (make-f64vector n 0.0)))
    (c_seulex_yy handle v n)
    v))


(define seulex-syy 
  (foreign-safe-lambda* scheme-object (((nonnull-c-pointer SEULEXSolverHandle) handle))
#<<EOF
    C_return (handle->syy);
EOF
))


(define (pointer+-f64 p n)
  (pointer+ p (fx* 8 n)))


)
