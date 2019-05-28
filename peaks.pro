pro	peaks,x,loci,y,fwhm=fwhm,threshold=threshold,plot=plot

;
;	Locates the position of peaks in an array x
;
;	IN: x 	  - float array		Array
; 	OUT: loci - float array		Positions for peaks in array x
;	        y - float array		Values in x at peaks
;
;	KEYWORDS: threshold - float	minimum value to be considered for
;					an acceptable peak
;		       plot - switch    when on it produces a plot to illustrate
;				        the peaks that have been found
;		       fwhm - float	expected fwhm of features (peaks); default is 3
;

if n_params() lt 2 then begin
	print,'% PEAKS: use -- peaks,x,loci,y[,fwhm=fwhm,threshold=threshold,plot=plot]' 
	return
endif

;default parameters
if n_elements(fwhm) eq 0 then fwhm=3. ; spected width of peaks
if n_elements(threshold) eq 0 then threshold=min(x)

;reset loci and y if need be
if n_elements(loci) gt 0 then undefine,loci
if n_elements(y) gt 0 then undefine,y

k=0
nx=n_elements(x)
if keyword_set(plot) then plot,x,yr=[0,max(x)*1.1],/ystyl
for i=1,nx-2 do begin
	if x[i] ge threshold and $
	x[i] ge max([x[max([0,fix(i-fwhm)]):i-1],x[i:min([nx-1,fix(i+fwhm)])]]) then begin
		c=poly_fit(findgen(3)-1.0,x[i-1:i+1],2)
		pox=float(i)-c[1]/2.d0/c[2]
		val=poly(-c[1]/2.d0/c[2],c)
		;print,i,pox
		;oplot,[i,i],[-10,10]
		if keyword_set(plot) then oplot,[pox,pox],[max(x)*1.02,max(x)*1.05],linestyle=2,thick=2
		if n_elements(loci) eq 0 then begin
			loci=pox
			y=val 
		endif else begin
			loci=[loci,pox]
			y=[y,val]
		endelse
	endif	
endfor

if keyword_set(plot) and n_elements(loci) gt 0 then xyouts,0.7,0.7,$
	'found '+string(n_elements(loci),format='(i4)')+' features',/norm

end

