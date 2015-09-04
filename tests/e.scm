;;
;; Compute Euler's number
;;

(use seulex srfi-4)

;; Problem Constants 

(define NEQ 1)

(define TEND  1.0)


(define (rhs t yy)
  (let ((dy (f64vector-ref yy 0)))
    (f64vector dy)))


(define (main)
  
  (let ((yy (f64vector 1.0))

	 ;; Integration limits 
	 (t0  0.0)
	 (tf  TEND)
	 (dt  1e-2))
    
    ;; solver initialization 
    (let ((solver (seulex-create-solver t0 yy rhs  
				     abstol: 1e-14
				     reltol: 1e-14)))

      ;; In loop, call solver, print results, and test for error. 
      
      (let recur ((tnext (+ t0 dt)) (iout 1))

	(let ((flag  (seulex-solve solver tnext)))
	  (if (negative? flag) (error 'main "SEULEX solver error" flag))

	  (print-results solver tnext)

	  (if (< tnext tf)
	      (recur (+ tnext dt) (+ 1 iout)))
	  ))

      (let ((yy (seulex-yy solver)))
	(assert (< (abs (- 2.71828182846 (f64vector-ref yy 0) )) 1e-12)) )

      
      (seulex-destroy-solver solver)
      
      )))



(define (print-results solver t)
  (let ((yy (seulex-yy solver)))
    (printf "~A ~A~%" 
	    t 
	    (f64vector-ref yy 0)
	    )))
      
      
(main)

