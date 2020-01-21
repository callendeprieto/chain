pro collapse1,f,idisp,ap1,s,vf=vf,vs=vs,width=width

;
; Collapse the spectrum assuming the dispersion is perfectly aligned with the detector
;
;	IN:      f - float array		Input image
;		idisp - integer (1/2)		Dispersion direction
;	        ap1 - float array		Locations of centers of the orders
;						in the cross-dispersion direction
;	OUT	 s - float array		Collapsed spectrum 
;	Keywords: fwhm	- float			approx. FWHM of the orders 
;		  width - float			ratio between the width of the spectral 
;						orders and the distance between them 
;						(0.72 is a good choice for HORS)
;		  smoothi
;		  smoothinglength -integer	length over which one should smooth the 
;						flux to help the peak search
;

sf=size(f)
np=sf[idisp] ; npixels in the dispersion direction
xdisp=1      
if idisp eq 1 then xdisp=2
nx=sf[xdisp] ; npixels in the cross-disp. direction

nap1=n_elements(ap1)

if not keyword_set(width) then width=0.72
delta=median(ap1-shift(ap1,1))/2.*width
s=fltarr(nap1,np)
for i=0,nap1-1 do begin
  s[i,*]=total(f[fix(ap1[i]-delta):fix(ap1[i]+delta),*],1)
endfor

if n_elements(vf) gt 0 then begin
	vs=fltarr(nap1,np)
	if idisp eq 2 then vf2=vf else vf2=transpose(vf) 
	for i=0,nap1-1 do begin
                vs[i,*]=total(vf2[fix(ap1[i]-delta):fix(ap1[i]+delta),*],1)
	endfor
endif

end
