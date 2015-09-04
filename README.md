# seulex

## Description

The Chicken `seulex` library provides bindings to the Fortran
`SEULEX` solver for systems of stiff differential and differential
algebraic equations of the form `MY'=F(X,Y)`. SEULEX was written by
E. Hairer AND G. Wanner and stands for Semi-implicit EULer
EXtrapolation and computes its error estimates by means of polynomial
extrapolation. Extrapolation solvers can achieve good performance for
certain types of problems where high precision is required.

## Installation notes

The Chicken `seulex` library must be compiled along with the
`SEULEX` Fortran source code. Therefore, the user must have a
Fortran compiler such as gfortran or g95 installed on their system.
The following environment variables can be used to control the Fortran
compilation process:

- `FORTRAN_COMPILER` : path to Fortran compiler (default is "gfortran")
- `FORTRAN_FLAGS` : flags to be passed to the Fortran compiler (default is "-fno-automatic -fPIC -g")
- `FORTRAN_LIBS` : Fortran libraries to link to  (default is "-lgfortran -lblas -llapack")


## Library procedures

<procedure>(seulex-create-solver XINIT YINIT FCN [JACOBIAN] [AUTONOMOUS] [USER-DATA] [DENSITY-COMPONENTS] [RELTOL] [ABSTOL]) => SEULEX-SOLVER</procedure>

Creates and initializes an object representing a problem to be solved
with the SEULEX solver.  

Arguments `XINIT` and `YINIT` represent the initial conditions:
`XINIT` is a scalar value for the independent variable, and
`YINIT` is an SRFI-4 `f64vector` value containing the initial
dependent variable values.

Argument `FCN` is used to compute the right-hand side function
`F` and must be a procedure of the following form:

 (LAMBDA T YY DATA)

or

 (LAMBDA T YY)

depending on whether the `USER-DATA` optional argument is set, where 

; `T`  :  real-valued independent variable
; `YY` : SRFI-4 `f64vector` with current variable values
; `DATA` : is a user data object (if set)

This procedure must return a SRFI-4 `f64vector` containing the
derivative values.

Optional keyword argument `JACOBIAN` is a procedure which will be
used to compute the partial derivatives of `F(X,Y)` with respect to
`Y`. If not given, it is computed internally by finite differences.

Optional keyword argument `AUTONOMOUS` is a boolean value that
indicates whether `F(X,Y)` is independent of `X` (autonomous) or
not (non-autonomous). The default is `#t` (autonomous).

Optional keyword argument `USER-DATA` is an object that will be
passed as an additional argument to the right-hand side function.

Optional keyword argument `DENSITY-COMPONENTS` must be an SRFI-4
`s32vector` which indicates the components for which dense output
is required.

Optional keyword arguments `RELTOL` and `ABSTOL` specify relative
and absolute error tolerance, respectively. These both default to
1e-4.

<procedure>(seulex-destroy-solver SEULEX-SOLVER)</procedure>

Deallocates the memory associated with the given solver.

<procedure>(seulex-solve SEULEX-SOLVER T)</procedure>

Integrates the system over an interval in the independent
variable. This procedure returns either when the given `T` is
reached, or when a root is found.

<procedure>(seulex-yy SEULEX-SOLVER)</procedure>

Returns the SRFI-4 `f64vector` value of current state values of the
system.

<procedure>(seulex-nfcn SEULEX-SOLVER)</procedure>

Returns the number of function evaluations done so far.

<procedure>(seulex-njac SEULEX-SOLVER)</procedure>

Returns the number of Jacobian function evaluations done so far.

<procedure>(seulex-nstep SEULEX-SOLVER)</procedure>

Returns the number of computed steps.

<procedure>(seulex-ndec SEULEX-SOLVER)</procedure>

Returns the number of LU decompositions.

<procedure>(seulex-nsol SEULEX-SOLVER)</procedure>

Returns the number of backward-forward substitutions.



## Example
>
> ;;
> ;; Hodgkin-Huxley model
> ;;
> 
> (use mathh seulex srfi-4)
> 
> (define neg -)
> (define pow expt)
> 
> (define TEND  500.0)
> 
>   	                   
> ;; Model parameters
> 
> (define (I_stim t) 10)
> (define C_m       1)
> (define E_Na      50)
> (define E_K       -77)
> (define E_L       -54.4)
>  (define gbar_Na   120)
> (define gbar_K    36)
> (define g_L       0.3)
> 
> ;; Rate functions
> 
> (define (amf v)   (* 0.1    (/ (+ v 40)  (- 1.0 (exp (/ (neg (+ v 40)) 10))))))
> (define (bmf v)   (* 4.0    (exp (/ (neg (+ v 65)) 18))))
> (define (ahf v)   (* 0.07   (exp (/ (neg (+ v 65)) 20))))
> (define (bhf v)   (/ 1.0    (+ 1.0 (exp (/ (neg (+ v 35)) 10)))))
> (define (anf v)   (* 0.01   (/ (+ v 55) (- 1 (exp (/ (neg (+ v 55)) 10))))))
> (define (bnf v)   (* 0.125  (exp (/ (neg (+ v 65)) 80))))
> 
> ;; State functions
> 
> (define (minf v) (* 0.5 (+ 1 (tanh (/ (- v v1) v2)))))
> (define (winf v) (* 0.5 (+ 1 (tanh (/ (- v v3) v4)))))
> (define (lamw v) (* phi (cosh (/ (- v v3) (* 2 v4)))))
>   	                   
> ;; Model equations
> 
> (define (rhs t yy)
> 
>   (let ((v (f64vector-ref yy 0))
> 	(m (f64vector-ref yy 1))
> 	(h (f64vector-ref yy 2))
> 	(n (f64vector-ref yy 3)))
> 
>     ;; transition rates at current step
>     (let ((am  (amf v))
> 	  (an  (anf v))
> 	  (ah  (ahf v))
> 	  (bm  (bmf v))
> 	  (bn  (bnf v))
> 	  (bh  (bhf v))
> 
> 	  (g_Na (* gbar_Na  (* h (pow m 3))))
> 	  (g_K  (* gbar_K   (pow n 4))))
>       
>       (let (
> 
> 	    ;; currents
> 	    (I_Na   (* (- v E_Na) g_Na))
> 	    (I_K    (* (- v E_K)  g_K))
> 	    (I_L    (* g_L  (- v E_L))))
> 		  
> 	(let (
> 	      ;; state equations
> 	      (dm (- (* am (- 1 m))  (* bm m)))
> 	      (dh (- (* ah (- 1 h))  (* bh h)))
> 	      (dn (- (* an (- 1 n))  (* bn n)))
> 	      (dv (/ (- (I_stim t) I_L I_Na I_K) C_m))
> 	      )
>    
> 	  (f64vector dv dm dh dn)
> 	  
> 	  )))
>     ))
>   
>  (let ((yy (f64vector -65  0.052 0.596 0.317)) ;; v m h n
> 
> 	;; Integration limits 
> 	(t0  0.0)
> 	(tf  TEND)
> 	(dt  1e-2))
>    
>     ;; solver initialization 
>     (let ((solver (seulex-create-solver
> 		   t0 yy rhs  
> 		   abstol: 1e-4
> 		   reltol: 1e-4)))
> 
>       ;; In loop, call solver, print results, and test for error. 
>       
>       (let recur ((tnext (+ t0 dt)) (iout 1))
> 
> 	(let ((flag  (seulex-solve solver tnext)))
> 	  (if (negative? flag) (error 'main "SEULEX solver error" flag))
>  
>          (print-results solver tnext)
> 
> 	  (if (< tnext tf)
> 	      (recur (+ tnext dt) (+ 1 iout)))
> 	  ))
>       
> 
>      (seulex-destroy-solver solver)
>       
> (define (print-results solver t)
>   (let ((yy (seulex-yy solver)))
>     (printf "~A ~A ~A ~A ~A~%" 
> 	    t 
> 	     (f64vector-ref yy 0)
>	     (f64vector-ref yy 1)
>	     (f64vector-ref yy 2)
>	     (f64vector-ref yy 3)
>	     )))
>

## Version History

- 1.0 Initial release

## License

>
>  The SEULEX solver was written by E. HAIRER AND G. WANNER, 
>  UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
> 
>  Chicken Scheme bindings for SEULEX are copyright 2011-2015 Ivan Raikov.
> 
> This program is free software: you can redistribute it and/or modify
> it under the terms of the GNU General Public License as published by
> the Free Software Foundation, either version 3 of the License, or (at
> your option) any later version.
> 
> This program is distributed in the hope that it will be useful, but
> WITHOUT ANY WARRANTY; without even the implied warranty of
> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
> General Public License for more details.
> 
> A full copy of the GPL license can be found at
> <http://www.gnu.org/licenses/>.
>


  
  
