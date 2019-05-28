pro	cut,f,datasec=datasec

;
;	Cut the image to remove overscan (or undesired sections of an image). 
;
;	IN/OUT:      f - float array		Input image
;       Keywords:   datasec - integer array      array with 4 values giving 
;						the coordinates of the corner
;						of the image to retain
;


sf=size(f)
f=float(f) ; in case it is not a float

;if n_elements(datasec) eq 0 then datasec=[1,4108,1,4094] ; HORS CCD
if n_elements(datasec) eq 0 then datasec=[1,4108,1,4088] ; HORS CCD
f=f[datasec[0]:datasec[1],datasec[2]:datasec[3]]

end

