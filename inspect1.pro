pro inspect1,f,idisp,left1,right1,x0=x0,dx=dx,y0=y0,dy=dy,_extra=e

;
; displays a 2d spectrum and a set of apertures assigned to it
;
;	IN:      f - float array		Input image
;		idisp - integer (1/2)		Dispersion direction
;	        left1   - int array		Array with location of left edges
;			(nap)	        for apertures
;	        right1   - int array		Array with location of right edges
;			(nap)	        for apertures
;	OUT	 none, just a plot
;      
;       KEYWORDS  
;                x0      int                    first pixel to show in x
;                                               (x0=0 by default)
;                dx      int                    range of pixels to show
;                                               (dx equals the full x range of
;                                               the image, if possible at all) 
;                y0      int                    first pixel to show in y
;                dy      int                    range of pixels to show in y
;                                               (dy equals the full y range of 
;                                               the image by default) 
;

sf=size(f)
np=sf[idisp] ; npixels in the dispersion direction
xdisp=1      
if idisp eq 1 then xdisp=2
nx=sf[xdisp] ; npixels in the cross-disp. direction

nap=n_elements(left1[*,0])

if n_elements(x0) eq 0 then x0=0
if n_elements(dx) eq 0 then dx=nx
if n_elements(y0) eq 0 then y0=0
if n_elements(dy) eq 0 then dy=np

display,f[x0:x0+dx-1,y0:y0+dy-1],/log,_extra=e
for i=0,nap-1 do begin
  oplot,replicate(left1[i],np)-x0,findgen(np)-y0,thick=2
  oplot,replicate(right1[i],np)-x0,findgen(np)-y0,thick=2,linestyle=2
endfor

end
