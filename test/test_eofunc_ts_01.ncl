begin
	nTime = 127;
	nLat  =  91;
	nLon  = 180;
	data = reshape(asciiread("resources/SampleSST.csv", -1, "double") , (/ nTime, nLat, nLon/))
	data!0 = "time"
	data!1 = "lat"
	data!2 = "lon"

	data_rearranged = data(lat|:, lon|:, time|:)

	neval = 5

	ev = eofunc(data_rearranged, neval, False)

	ev_ts = eofunc_ts(data_rearranged, ev, False)

	printVarSummary(ev)
	printVarSummary(ev_ts)

	asciiwrite("resources/SampleSST_ev_rearranged.txt", ev)
	asciiwrite("resources/SampleSST_ev_rearranged_ts.txt", ev_ts)
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

;Variable: ev_ts
;Type: double
;Total Size: 5080 bytes
;            635 values
;Number of Dimensions: 2
;Dimensions and sizes:	[5] x [127]
;Coordinates: 
;Number Of Attributes: 3
;  _FillValue :	9.969209968386869e+36
;  matrix :	covariance
;  ts_mean :	( 258.2679041145958, 381.9555081707763, 257.2192800344826, 611.9719389444721, -193.1237981848278 )