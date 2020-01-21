pro inspect1,f,idisp,ap1,width=width,_extra=e

;
; displays a 2d spectrum and a set of apertures assigned to it
;
;	IN:      f - float array		Input image
;		idisp - integer (1/2)		Dispersion direction
;	        ap1 - float array		Locations of centers of the orders
;						in the cross-dispersion direction
;	OUT	 none, just a plot
;
;	Keywords:
;		  width - float			ratio between the width of the spectral 
;						orders and the distance between them 
;						(0.72 is a good choice for HORS)
;

sf=size(f)
np=sf[idisp] ; npixels in the dispersion direction
xdisp=1      
if idisp eq 1 then xdisp=2
nx=sf[xdisp] ; npixels in the cross-disp. direction

nap1=n_elements(ap1)

if not keyword_set(width) then width=0.72
delta=median(ap1-shift(ap1,1))/2.*width
print,'delta(inspect1)=',delta
display,f,/log,_extra=e
oplot,ap1,replicate(np/2,n_elements(ap1)),psy=2,thick=3
for i=0,n_elements(ap1)-1 do begin
  oplot,[fix(ap1[i]-delta),fix(ap1[i]-delta)],[0,1e5],thick=2
  oplot,[fix(ap1[i]+delta),fix(ap1[i]+delta)],[0,1e5],thick=2,linestyle=2
endfor

end
