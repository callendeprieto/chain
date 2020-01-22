pro collapse1,f,idisp,left1,right1,s,vf=vf,vs=vs

;
; Collapse the spectrum assuming the dispersion is perfectly aligned with the detector
;
;	IN:      f - float array		Input image
;		idisp - integer (1/2)		Dispersion direction
;	        left1   - float array		Array with location of left edges
;			(nap)	                for apertures
;	        right1   - float array		Array with location of right edges
;			(nap)	                for apertures

;	OUT	 s - float array		Collapsed spectrum 
;	Keywords: fwhm	- float			approx. FWHM of the orders 
;		  smoothi
;		  smoothinglength -integer	length over which one should smooth the 
;						flux to help the peak search
;

sf=size(f)
np=sf[idisp] ; npixels in the dispersion direction
xdisp=1      
if idisp eq 1 then xdisp=2
nx=sf[xdisp] ; npixels in the cross-disp. direction

nap1=n_elements(left1)

s=fltarr(nap1,np)
for i=0,nap1-1 do begin
  s[i,*]=total(f[left1[i]:right1[i],*],1)
endfor

if n_elements(vf) gt 0 then begin
	vs=fltarr(nap1,np)
	if idisp eq 2 then vf2=vf else vf2=transpose(vf) 
	for i=0,nap1-1 do begin
                vs[i,*]=total(vf2[left1[i]:right1[i],*],1)
	endfor
endif

end
