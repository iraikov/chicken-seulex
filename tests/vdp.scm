;;
;; Van Der Pol oscillator
;;

(use mathh seulex srfi-4)

(define neg -)
(define pow expt)

(define TEND  2.0)

  	                   
;; Model parameters

(define r       1e-6)

;; Right-hand side of equations

(define (rhs t yy)

  (let ((x (f64vector-ref yy 0))
	(y (f64vector-ref yy 1)))

    (let ((dx y)
	  (dy (/ (- (* (- 1 (* x x)) y) x) r)))
    
      (f64vector dx dy)))
  )

;; Jacobian (column-major order)

(define (vdpjac t yy)

  (let ((x (f64vector-ref yy 0))
	(y (f64vector-ref yy 1)))
    
    (let ((df11 0.0) (df12 1.0)
	  (df21 (/ (- (* -2.0 x y) 1.0) r))
	  (df22 (/ (- 1.0 (* x x)) r)))

      (f64vector df11 df21 df12 df22)

      )))


(define (main)
  
  (let ((yy (f64vector 2.0 -0.066)) 

	;; Integration limits 
	(t0  0.0)
	(tf  TEND)
	(dt  1e-2))
    
    ;; solver initialization 
    (let ((solver (seulex-create-solver
		   t0 yy rhs  
		   jacobian: vdpjac
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
    (printf "~A ~A ~A~%" 
	    t 
	    (f64vector-ref yy 0)
	    (f64vector-ref yy 1)
	    )))
      
      
(main)
