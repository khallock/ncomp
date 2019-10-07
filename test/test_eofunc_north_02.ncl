begin
	eval = (/ 974.8881, 490.3863, 190.7456, 172.2483 /)

	sig = eofunc_north(eval, 80, True)

	print(sig)
end

;(0)	154.143   820.745   974.888   1129.03  True
;(1)	77.5369   412.849   490.386   567.923  True
;(2)	30.1595   160.586   190.746   220.905  False
;(3)	27.2348   145.013   172.248   199.483  False


;Variable: sig
;Type: logical
;Total Size: 16 bytes
;            4 values
;Number of Dimensions: 1
;Dimensions and sizes:	[4]
;Coordinates: 
;Number Of Attributes: 3
;  _FillValue :	Missing
;  long_name :	EOF separation
;  N :	80
;(0)	True
;(1)	True
;(2)	False
;(3)	False