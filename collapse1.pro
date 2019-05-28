pro collapse1,f,idisp,ap1,s

;
; Collapse the spectrum assuming the dispersion is perfectly aligned with the detector
;
;	IN:      f - float array		Input image
;		idisp - integer (1/2)		Dispersion direction
;	        ap1 - float array		Locations of centers of the orders
;						in the cross-dispersion direction
;	OUT	 s - float array		Collapsed spectrum 
;	Keywords: fwhm	- float			approx. FWHM of the orders 
;		  smoothinglength -integer	length over which one should smooth the 
;						flux to help the peak search
;

sf=size(f)
np=sf[idisp] ; npixels in the dispersion direction
xdisp=1      
if idisp eq 1 then xdisp=2
nx=sf[xdisp] ; npixels in the cross-disp. direction

nap1=n_elements(ap1)

delta=median(ap1-shift(ap1,1))/2.*0.7
s=fltarr(nap1,np)
for i=0,nap1-1 do begin
  s[i,*]=total(f[fix(ap1[i]-delta):fix(ap1[i]+delta),*],1)
endfor

end
