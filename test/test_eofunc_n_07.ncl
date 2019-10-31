begin
	nTime = 127;
	nLat  =  91;
	nLon  = 180;
	data = reshape(asciiread("resources/SampleSST.csv", -1, "double") , (/ nTime, nLat, nLon/))

	neval = 5
	option = False

	ev = eofunc_n(data, neval, option, 0)

	printVarSummary(ev)

	asciiwrite("SampleSST_ev.txt", ev)
end


;Variable: ev
;Type: double
;Total Size: 655200 bytes
;            81900 values
;Number of Dimensions: 3
;Dimensions and sizes:	[5] x [91] x [180]
;Coordinates: 
;Number Of Attributes: 6
;  _FillValue :	9.969209968386869e+36
;  method :	transpose
;  matrix :	covariance
;  pcvar :	( 86.40701, 3.54013, 2.655234, 1.046697, 0.7760174 )
;  eval :	( 61037.92580653914, 2500.748587296903, 1875.65759356925, 739.3867320896884, 548.1788335801401 )
;  eval_transpose :	( 570.4479047340106, 23.37148212427012, 17.52951022027337, 6.910156374669985, 5.123166668973272 )