pro collapse,f,idisp,left,right,s,vs=vs,vf=vf

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
;

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
  for j=0,np-1 do vs[i,j]=total(vf2[left[i,j]:right[i,j],j])
endfor

for i=0,nap-1 do begin
  for j=0,np-1 do begin
    s[i,j]=total(f2[left[i,j]:right[i,j],j])
  endfor
endfor



end
