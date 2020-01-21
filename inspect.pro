pro inspect,f,idisp,left,right,_extra=e

;
; displays a 2d spectrum and a set of apertures assigned to it
;
;	IN:      f - float array		Input image
;		idisp - integer (1/2)		Dispersion direction
;	        left   - int array		Array with location of left edge
s
;			(nap x np)	        for apertures
;	        right   - int array		Array with location of right edg
es
;			(nap x np)	        for apertures
;	OUT	 none, just a plot
;

sf=size(f)
np=sf[idisp] ; npixels in the dispersion direction
xdisp=1      
if idisp eq 1 then xdisp=2
nx=sf[xdisp] ; npixels in the cross-disp. direction

nap=n_elements(ap[*,0])


display,f,/log,_extra=e
oplot,ap[*,np/2],replicate(np/2,nap),psy=2,thick=3
for i=0,nap-1 do begin
  oplot,left[i,*],findgen(np),thick=2
  oplot,right[i,*],findgen(np),thick=2,linestyle=2
endfor

end
