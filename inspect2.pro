pro inspect2,f,idisp,ap,width=width,_extra=e

;
; displays a 2d spectrum and a set of apertures assigned to it
;
;	IN:      f - float array		Input image
;		idisp - integer (1/2)		Dispersion direction
;                ap   - float array		Array with location of flux-weighted
;			(nap x np)	        aperture centers
;	OUT	 none, just a plot
;	Keywords:
;		  width - float			ratio between the width of the spectral 
;						orders and the distance between them 
;						(0.72 is a good choice for HORS)

sf=size(f)
np=sf[idisp] ; npixels in the dispersion direction
xdisp=1      
if idisp eq 1 then xdisp=2
nx=sf[xdisp] ; npixels in the cross-disp. direction

nap=n_elements(ap[*,0])

if not keyword_set(width) then width=0.72
delta=median(ap[*,np/2]-shift(ap[*,np/2],1))/2.*width
print,'delta(inspect2)=',delta


display,f,/log,_extra=e
oplot,ap[*,np/2],replicate(np/2,nap),psy=2,thick=3
for i=0,nap-1 do begin
  oplot,fix(ap[i,*]+delta)-fix(2.*delta)-1,findgen(np),thick=2
  oplot,fix(ap[i,*]+delta+1),findgen(np),thick=2,linestyle=2
endfor

end
