pro inspect,f,idisp,ap,delta1

;
; displays a 2d spectrum and a set of apertures assigned to it
;
;	IN:      f - float array		Input image
;		idisp - integer (1/2)		Dispersion direction
;                ap   - float array		Array with location of flux-weighted
;			(nap x np)	        aperture centers
;	delta1     - float                      Half width of the orders in the xdisp dir
;						(extraction should proceed between the
;						centers in the ap array and -/+ delta1)
;	OUT	 none, just a plot

sf=size(f)
np=sf[idisp] ; npixels in the dispersion direction
xdisp=1      
if idisp eq 1 then xdisp=2
nx=sf[xdisp] ; npixels in the cross-disp. direction

nap=n_elements(ap[*,0])

display,f,/log
oplot,ap[*,np/2],replicate(np/2,nap),psy=2
for i=0,nap-1 do begin
  oplot,ap[i,*]-delta1,findgen(np),thick=2,linestyle=1
  oplot,ap[i,*]+delta1,findgen(np),thick=2,linestyle=0
endfor

end
