pro collapse,f,idisp,left,right,s,vs=vs,vf=vf,dleft=dleft,dright=dright,clean=clean

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
;		 vs   - float array		Variance spectrum (when vf is not made available, 
; 						 it is approximated vf=f)
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

;print,'p1'
;help,/mem


if keyword_set(clean) then begin

  ;spatial profile for weighting
  mm=minmax(right-left+1)
  if mm[0] ne mm[1] then begin
    print,'% COLLAPSE: the width of the extraction window is not constant, cannot clean!'
    stop
  endif else begin
    p=fltarr(nap,mm[0]) ; spatial profile
    ;a=fltarr(nap)         ; collapsed profile
  endelse

  ;print,'p2'
  ;help,/mem


  for i=0,nap-1 do begin
    for j=0,np-1 do begin
      p[i,*]=p[i,*]+f2[left[i,j]+dleft:right[i,j]+dright,j]
    endfor
    ;a[i]=total(p[i,*])
  endfor
  p=p/np
  ;a=a/np
  ;w=1./p

  ;print,'p3'
  help,/mem

  for i=0,nap-1 do begin
    for j=0,np-1 do begin
      arr=f2[left[i,j]+dleft:right[i,j]+dright,j]
      me=median(arr/p[i,*]) ; median factor (profile at this wave)/(average prof)
      ma=mad(arr/p[i,*]-me) ; mad
      wbad=where(arr/p[i,*]-me gt 20.*ma)
      if max(wbad) gt -1 then arr[wbad]=p[i,wbad]*me
      s[i,j]=total(arr) 
      ; total(arr*w[i,*])/total(w[i,*])*mm[0]
      arr=vf2[left[i,j]+dleft:right[i,j]+dright,j]
      vs[i,j]=total(arr) 
      ; total(vf2[left[i,j]+dleft:right[i,j]+dright,j]*w[i,*])/total(w[i,*])*mm[0]
    endfor
  endfor
  
  ;print,'p4'
  help,/mem


endif else begin

  for i=0,nap-1 do begin
    for j=0,np-1 do begin
      s[i,j]=total(f2[left[i,j]+dleft:right[i,j]+dright,j])
      vs[i,j]=total(vf2[left[i,j]+dleft:right[i,j]+dright,j])
    endfor
  endfor

endelse

;stop

end
