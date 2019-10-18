begin
	m = todouble( (/ 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63 /) )

	n = reshape(m, (/ 4,4,4 /))

	opt = False
	opt@debug = True

	r = eofunc_n(n, 1, opt, 1)

	print(r)
end


;Variable: r
;Type: double
;Total Size: 128 bytes
;            16 values
;Number of Dimensions: 3
;Dimensions and sizes:	[1] x [4] x [4]
;Coordinates: 
;Number Of Attributes: 6
;  _FillValue :	9.969209968386869e+36
;  method :	transpose
;  matrix :	covariance
;  pcvar :	100
;  eval :	426.6666666666666
;  eval_transpose :	85.33333333333333
;(0,0,0)	0.25
;(0,0,1)	0.25
;(0,0,2)	0.25
;(0,0,3)	0.25
;(0,1,0)	0.25
;(0,1,1)	0.25
;(0,1,2)	0.25
;(0,1,3)	0.25
;(0,2,0)	0.25
;(0,2,1)	0.25
;(0,2,2)	0.25
;(0,2,3)	0.25
;(0,3,0)	0.25
;(0,3,1)	0.25
;(0,3,2)	0.25
;(0,3,3)	0.25