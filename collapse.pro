pro collapse,f,idisp,left,right,s,vs=vs,vf=vf,dleft=dleft,dright=dright

;
; Collapse the spectrum
;
;	IN:      f - float array		Input image
;		idisp - integer (1/2)		Dispersion direction
;	        left   - float array		Array with location of left edges
;			(nap x np)	        for apertures
;	        right   - float array		Array with location of right edges
;			(nap x np)	        for apertures
;
;	OUT	 s - float array		Collapsed spectrum 
;	KEYWORD: vf   - float array		Variance image (image)
;		 vs   - float array		Variance spectrum(when vf is made available)
;                dleft - int                    shift left aperture limits by these many pixels
;                dright - int                   shift right aperture limits by these many pixels
;                                               (default is dleft=dright=0, but this allows extraction
;                                               of individual fibers (e.g. use dleft=0 & dright=-10 to extract
;                                               the first fiber/pixel on the left of the apertures)
;

if n_elements(dleft) eq 0 then dleft=0
if n_elements(dright) eq 0 then dright=0
sf=size(f)
np=sf[idisp] ; npixels in the dispersion direction
xdisp=1      
if idisp eq 1 then xdisp=2
nx=sf[xdisp] ; npixels in the cross-disp. direction

nap=n_elements(left[*,0])
s=fltarr(nap,np)
if idisp eq 2 then f2=f else f2=transpose(f)


if n_elements(vf) eq 0 then vf=f
vs=fltarr(nap,np)
if idisp eq 2 then vf2=vf else vf2=transpose(vf)

for i=0,nap-1 do begin
  for j=0,np-1 do vs[i,j]=total(vf2[left[i,j]+dleft:right[i,j]+dright,j])
endfor

for i=0,nap-1 do begin
  for j=0,np-1 do begin
    s[i,j]=total(f2[left[i,j]+dleft:right[i,j]+dright,j])
  endfor
endfor



end
