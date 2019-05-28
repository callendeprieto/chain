pro collapse2,f,idisp,ap,delta1,s,vs=vs,vf=vf

;
; Collapse the spectrum
;
;	IN:      f - float array		Input image
;		idisp - integer (1/2)		Dispersion direction
;	        ap   - float array		Array with location of flux-weighted
;			(nap x np)	        aperture centers
;	delta1     - float                      Half width of the orders in the xdisp dir
;						(extraction should proceed between the
;						centers in the ap array and -/+ delta1)
;
;	OUT	 s - float array		Collapsed spectrum 
;	KEYWORD: vf   - float array		Variance image (image)
;		 vs   - float array		Variance spectrum(when vf is made available)
;
;

sf=size(f)
np=sf[idisp] ; npixels in the dispersion direction
xdisp=1      
if idisp eq 1 then xdisp=2
nx=sf[xdisp] ; npixels in the cross-disp. direction

nap=n_elements(ap[*,0])
s=fltarr(nap,np)
if idisp eq 2 then f2=f else f2=transpose(f)


if n_elements(vf) eq 0 then vf=f
vs=fltarr(nap,np)
if idisp eq 2 then vf2=vf else vf2=transpose(vf)
for i=0,nap-1 do begin
  for j=0,np-1 do vs[i,j]=total(vf2[fix(ap[i,j]-delta1):fix(ap[i,j]+delta1),j])
endfor

for i=0,nap-1 do begin
  ;for j=0,np-1 do s[i,j]=total(f2[fix(ap[i,j]-delta1):fix(ap[i,j]+delta1),j])
  for j=0,np-1 do begin
	data=f2[fix(ap[i,j]-delta1):fix(ap[i,j]+delta1),j]
	vdata=vf2[fix(ap[i,j]-delta1):fix(ap[i,j]+delta1),j]
	frac=data/total(data)
	den=total(frac^2/vdata)
	w=frac/vdata/den
	s[i,j]=total(w*f2[fix(ap[i,j]-delta1):fix(ap[i,j]+delta1),j])
  endfor
endfor



end
