;;
;; Morris-Lecar model
;;

(use mathh seulex srfi-4)


(define TEND  5.0)

  	                   
;; Model parameters

(define Istim  50.0)
(define vl     -50)
(define vk     -70)
(define vca    100)
(define gl     2.0)
(define gk     8.0)
(define gca    4.0)
(define c      20.0)
(define v1     -1.0)
(define v2     15)
(define v3     10)
(define v4     14.5)
(define phi    0.0667)


;; State functions

(define (minf v) (* 0.5 (+ 1 (tanh (/ (- v v1) v2)))))
(define (winf v) (* 0.5 (+ 1 (tanh (/ (- v v3) v4)))))
(define (lamw v) (* phi (cosh (/ (- v v3) (* 2 v4)))))
  	                   
;; Model equations

(define (rhs t yy)


  (let ((v (f64vector-ref yy 0))
	(w (f64vector-ref yy 1)))

  (let ((ica (* gca (minf v)  (- vca v)))
	(ik  (* gk w (- vk v ))))
    
    (let ((dv (/ (+ Istim (* gl (- vl v)) ica ik) c))
	  (dw (* (lamw v) (- (winf v) w))))

      (f64vector dv dw)

      ))
  ))


(define (main)
  
  (let ((yy (f64vector -60.899 0.0149));; v w

	;; Integration limits 
	(t0  0.0)
	(tf  TEND)
	(dt  1))
    
    ;; solver initialization 
    (let ((solver (seulex-create-solver
		   t0 yy rhs
		   abstol: 1e-4
		   reltol: 1e-4)))

      ;; In loop, call the solver, print results, and test for error. 

      (let recur ((tnext (+ t0 dt))
		  (iout 1))

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
