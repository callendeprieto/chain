pro	hbias,f,overscan=overscan,rdn=rdn,bin=bin

;
;	Bias removal using the overscan area. We subtract the median 
;	signal. 
;
;	IN/OUT:      f - float array		Input image
;	Keywords:  overscan - integer array 	array with 3 values
;						1) readout direction (1/2)
;						2,3) overscan start/end pixels  
;		     rdn  - float		estimated readout noise
;
;		     bin - integer		set to indicate 8x2 binning
;
;


sf=size(f)

f=float(f) ; in case it is not a float

if n_elements(overscan) eq 0 and keyword_set(bin) then overscan=[1,1,1]
if n_elements(overscan) eq 0 then overscan=[1,1,16] ; HORuS CCD

if overscan[0] eq 1 then begin
	f=f-median(f[overscan[1]:overscan[2],*])
endif else begin
	f=f-median(f[*,overscan[1]:overscan[2]])
endelse

if overscan[0] eq 1 then rdn=stddev(f[overscan[1]:overscan[2],*]) else $
	rdn=stddev(f[*,overscan[1]:overscan[2]])


end

