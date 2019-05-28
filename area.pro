pro area,f,idisp,ap,delta1,w

;
; return indices of the orders
;
;

sf=size(f)
np=sf[idisp] ; npixels in the dispersion direction
xdisp=1      
if idisp eq 1 then xdisp=2
nx=sf[xdisp] ; npixels in the cross-disp. direction

w=bytarr(nx,np)
nap=n_elements(ap[*,0])
s=fltarr(nap,np)
if idisp eq 2 then f2=f else f2=transpose(f) 
for i=0,nap-1 do begin
  for j=0,np-1 do w[fix(ap[i,j]-delta1):fix(ap[i,j]+delta1),j]=1
endfor



end
