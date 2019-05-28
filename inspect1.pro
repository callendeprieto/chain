pro inspect1,f,idisp,ap1

;
; displays a 2d spectrum and a set of apertures assigned to it
;
;	IN:      f - float array		Input image
;		idisp - integer (1/2)		Dispersion direction
;	        ap1 - float array		Locations of centers of the orders
;						in the cross-dispersion direction
;	OUT	 none, just a plot

sf=size(f)
np=sf[idisp] ; npixels in the dispersion direction
xdisp=1      
if idisp eq 1 then xdisp=2
nx=sf[xdisp] ; npixels in the cross-disp. direction

nap1=n_elements(ap1)

delta=median(ap1-shift(ap1,1))/2.*0.7
display,f,/log
oplot,ap1,replicate(np/2,n_elements(ap1)),psy=2
for i=0,n_elements(ap1)-1 do begin
  oplot,[ap1[i]-delta,ap1[i]-delta],[0,1e5]
  oplot,[ap1[i]+delta,ap1[i]+delta],[0,1e5]
endfor

end
