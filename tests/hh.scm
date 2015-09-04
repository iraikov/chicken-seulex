;;
;; Hodgkin-Huxley model
;;

(use mathh seulex srfi-4)

(define neg -)
(define pow expt)

(define TEND  5.0)

  	                   
;; Model parameters

;; model parameters
(define (I_stim t) 10)
(define C_m       1)
(define E_Na      50)
(define E_K       -77)
(define E_L       -54.4)
(define gbar_Na   120)
(define gbar_K    36)
(define g_L       0.3)

;; Rate functions

(define (amf v)   (* 0.1    (/ (+ v 40)  (- 1.0 (exp (/ (neg (+ v 40)) 10))))))
(define (bmf v)   (* 4.0    (exp (/ (neg (+ v 65)) 18))))
(define (ahf v)   (* 0.07   (exp (/ (neg (+ v 65)) 20))))
(define (bhf v)   (/ 1.0    (+ 1.0 (exp (/ (neg (+ v 35)) 10)))))
(define (anf v)   (* 0.01   (/ (+ v 55) (- 1 (exp (/ (neg (+ v 55)) 10))))))
(define (bnf v)   (* 0.125  (exp (/ (neg (+ v 65)) 80))))

;; State functions

(define (minf v) (* 0.5 (+ 1 (tanh (/ (- v v1) v2)))))
(define (winf v) (* 0.5 (+ 1 (tanh (/ (- v v3) v4)))))
(define (lamw v) (* phi (cosh (/ (- v v3) (* 2 v4)))))
  	                   
;; Model equations

(define (rhs t yy)

  (let ((v (f64vector-ref yy 0))
	(m (f64vector-ref yy 1))
	(h (f64vector-ref yy 2))
	(n (f64vector-ref yy 3)))

    ;; transition rates at current step
    (let ((am  (amf v))
	  (an  (anf v))
	  (ah  (ahf v))
	  (bm  (bmf v))
	  (bn  (bnf v))
	  (bh  (bhf v))

	  (g_Na (* gbar_Na  (* h (pow m 3))))
	  (g_K  (* gbar_K   (pow n 4))))
      
      (let (

	    ;; currents
	    (I_Na   (* (- v E_Na) g_Na))
	    (I_K    (* (- v E_K)  g_K))
	    (I_L    (* g_L  (- v E_L))))
		  
	(let (
	      ;; state equations
	      (dm (- (* am (- 1 m))  (* bm m)))
	      (dh (- (* ah (- 1 h))  (* bh h)))
	      (dn (- (* an (- 1 n))  (* bn n)))
	      (dv (/ (- (I_stim t) I_L I_Na I_K) C_m))
	      )
    
	  (f64vector dv dm dh dn)
	  
	  )))
    ))

(define (main)
  
  (let ((yy (f64vector -65  0.052 0.596 0.317)) ;; v m h n

	;; Integration limits 
	(t0  0.0)
	(tf  TEND)
	(dt  1e-2))
    
    ;; solver initialization 
    (let ((solver (seulex-create-solver
		   t0 yy rhs  
		   abstol: 1e-4
		   reltol: 1e-4)))

      ;; In loop, call the solver, print results, and test for error. 
      
      (let recur ((tnext (+ t0 dt)) (iout 1))

	(let ((flag  (seulex-solve solver tnext)))
	  (if (negative? flag) (error 'main "SEULEX solver error" flag))

	  (print-results solver tnext)

	  (if (< tnext tf)
	      (recur (+ tnext dt) (+ 1 iout)))
	  ))
      

      
      (seulex-destroy-solver solver)
      
      )))


(define (print-results solver t)
  (let ((yy (seulex-yy solver)))
    (printf "~A ~A ~A ~A ~A~%" 
	    t 
	    (f64vector-ref yy 0)
	    (f64vector-ref yy 1)
	    (f64vector-ref yy 2)
	    (f64vector-ref yy 3)
	    )))
      
      
(main)
