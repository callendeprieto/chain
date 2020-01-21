pro collapse2,f,idisp,ap,s,vs=vs,vf=vf,width=width

;
; Collapse the spectrum
;
;	IN:      f - float array		Input image
;		idisp - integer (1/2)		Dispersion direction
;	        ap   - float array		Array with location of flux-weighted
;			(nap x np)	        aperture centers
;
;	OUT	 s - float array		Collapsed spectrum 
;	KEYWORD: vf   - float array		Variance image (image)
;		 vs   - float array		Variance spectrum(when vf is made available)
;		  width - float			ratio between the width of the spectral 
;						orders and the distance between them 
;						(0.72 is a good choice for HORS)
;
;

sf=size(f)
np=sf[idisp] ; npixels in the dispersion direction
xdisp=1      
if idisp eq 1 then xdisp=2
nx=sf[xdisp] ; npixels in the cross-disp. direction

nap=n_elements(ap[*,0])
s=fltarr(nap,np)
sk=s
if idisp eq 2 then f2=f else f2=transpose(f)

if not keyword_set(width) then width=0.72
delta=median(ap[*,np/2]-shift(ap[*,np/2],1))/2.*width

if n_elements(vf) eq 0 then vf=f
vs=fltarr(nap,np)
if idisp eq 2 then vf2=vf else vf2=transpose(vf)
for i=0,nap-1 do begin
  for j=0,np-1 do vs[i,j]=total(vf2[fix(ap[i,j]+delta)-fix(2*delta)-1:fix(ap[i,j]+delta+1),j])
endfor

for i=0,nap-1 do begin
  for j=0,np-1 do begin
    s[i,j]=total(f2[fix(ap[i,j]+delta)-fix(2*delta)-1:fix(ap[i,j]+delta+1),j])
  endfor
;  for j=0,np-1 do begin
;	data=f2[fix(ap[i,j]+delta)-fix(2*delta)-1:fix(ap[i,j]+delta),j]
;	vdata=vf2[fix(ap[i,j]+delta)-fix(2*delta)-1:fix(ap[i,j]+delta),j]
;	frac=data/total(data)
;	den=total(frac^2/vdata)
;	w=frac/vdata/den
;	s[i,j]=total(w*f2[fix(ap[i,j]+delta)-fix(2*delta)-1:fix(ap[i,j]+delta),j])
;  endfor
endfor



end
