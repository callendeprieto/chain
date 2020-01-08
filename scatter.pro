pro  scatter,f,idisp,delta1,back
;+
; Removal of scattered light between orders
;
; We determine, for each column in the detector in the direction perpendicular to 
; the dispersion. a running 8 percentile with a width of the same as the orders 
; full width (2*delta1), and smooth it with a running mean using a width 3 times 
; larger. After smoothing the result in the dispersion direction, we return the result
; in the 'back' image.
;
; IN: f -- 2D spectral image
;     idisp -- direction of dispersion (1 along columns, 2 along rows)
;     delta1 -- half-width of the orders in the spatial direction
; 
; OUT: back -- 2D image with the estimated scattered light background
;-

if N_Params() lt 4 then begin
        print,'% scatter: use -- scatter,im,idisp,delta1,back'
	return
endif

sf=size(f)
np=sf[idisp] ; npixels in the dispersion direction
xdisp=1      
if idisp eq 1 then xdisp=2
nx=sf[xdisp] ; npixels in the cross-disp. direction

if idisp eq 1 then f=transpose(f)
back=f
back[*,*]=0.0

for i=0,np-1 do begin
  res=cmedian(f[*,i],delta1*3.,perc=0.08)
  back[*,i]=csmooth(res,delta1*6)
endfor

for i=0,nx-1 do begin
  if idisp eq 2 then res=transpose(back[i,*]) else res=back[i,*]
  back[i,*]=csmooth(res,delta1)
endfor

if idisp eq 1 then begin
  f=transpose(f)
  back=transpose(back)
endif

end
