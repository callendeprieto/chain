pro traze,f,idisp,ap1,ap,$
width=width,smoothinglength=smoothinglength,order=order,ntrace=ntrace,first=first,last=last

;
; Trace orders in the dispersion direction
;
;	IN:      f - float array		Input image
;		idisp - integer (1/2)		Dispersion direction
;	        ap1 - float array		Locations of centers of the orders
;						in the cross-dispersion direction
;       OUT:     ap   - float array		Array with location of flux-weighted
;			(nap x np)	        aperture centers
;	Keywords:
;		  width - float			ratio between the width of the spectral 
;						orders and the distance between them 
;						(0.7 is a good choice for HORS)
;		  smoothinglength -integer	length over which one should smooth the 
;						flux to help the peak search
;		  order  - integer		order for the polynomial fitting to the 
;					        distortion of the spectral orders 
;						(deviations from the perfect
;						alignment of the orders with the detector)
;		  ntrace - integer		Number of points to trace order curvature
;
;		  first - integer		Set this to the number of the first aperture to 
;						consider in the tracing. 
;
;		  last - integer		Set this to the number of the last apertures to 
;						consider in the tracing. 
;
;						Before 'first' and after 'last' the orders will use the 
;						same polynomial derived from first and last, respectively, 
;						but shifted according to the aperture centers in ap1
;

sf=size(f)
np=sf[idisp] ; npixels in the dispersion direction
xdisp=1      
if idisp eq 1 then xdisp=2
nx=sf[xdisp] ; npixels in the cross-disp. direction


if not keyword_set(width) then width=0.7
if not keyword_set(smoothinglength) then smoothinglength=4
if n_elements(order) eq 0 then order=2
if not keyword_set(ntrace) then ntrace=40

nap1=n_elements(ap1)
if not keyword_set(first) then first=0
if not keyword_set(last) then last=nap1-1
ap=fltarr(nap1,np)
delta=median(ap1-shift(ap1,1))/2.*width

for i=first,last do begin

  ;if i gt 0 then delta=(ap1[i]-ap1[i-1])/2.*width else $
  ;	delta=(ap1[i+1]-ap1[i])/2.*width


  delta=median(ap1-shift(ap1,1))/2.*width

  tpos=fltarr(ntrace)
  tcent=fltarr(ntrace)
  ;search for offsets, from center up
  pcent=ap1[i]
  for j=ntrace/2,ntrace-1 do begin
   tpos[j]=np/ntrace*(2.*j+1.)/2.
   prof=total(f[fix(pcent-delta):fix(pcent+delta),j*np/ntrace:(j+1)*np/ntrace-1],2)
   ;prof=prof-min(prof)
   cent=total(findgen(n_elements(prof))*prof)/total(prof)
   ;cent=median(findgen(n_elements(prof))*prof)/median(prof)
   tcent[j]=pcent+cent-delta
   ;print,ap1[i],pcent,tcent[j]
   pcent=tcent[j]
   ;plots,tcent[j],tpos[j],psy=4
   ;print,'ap #',i,'piece #',j
   ;plot,prof
   ;oplot,[cent,cent],[0,max(prof)]   
   ;oplot,[ap1[i]-pcent+delta,ap1[i]-pcent+delta],[0,max(prof)],thick=2,linestyle=2

	;stop

   ;if i gt 20 then stop
  endfor

  ;from center down
  pcent=ap1[i]
  for j=ntrace/2-1,0,-1 do begin
   tpos[j]=np/ntrace*(2.*j+1.)/2.
   prof=total(f[fix(pcent-delta):fix(pcent+delta),j*np/ntrace:(j+1)*np/ntrace-1],2)
   ;cent=total(findgen(n_elements(prof))*prof)/total(prof)
   cent=median(findgen(n_elements(prof))*prof)/median(prof)
   tcent[j]=pcent+cent-delta
   ;print,ap1[i],pcent,tcent[j]
   pcent=tcent[j]
   ;plots,tcent[j],tpos[j],psy=4
   ;print,'ap #',i,'piece #',j
   ;plot,prof
   ;oplot,[cent,cent],[0,max(prof)]   
   ;oplot,[ap1[i]-pcent+delta,ap1[i]-pcent+delta],[0,max(prof)],thick=2,linestyle=2
   ;print,tpos
   ;print,tcent
   ;if i gt 20 then stop
  endfor

  ;plot,tpos,tcent,psy=-2,yr=[min(tcent),max(tcent)],/ystyl
  c=poly_fit(tpos,csmooth(tcent,smoothinglength),order)
  ;oplot,tpos,poly(tpos,c),thick=2,col=180
  ap[i,*]=poly(findgen(np),c)
  ;print,c
  ;help,order,c
  ;wait,0.5
  ;ap[i,*]=interpol(tcent,tpos,findgen(np))

endfor

if first gt 0 then begin
	for i=0,first-1 do begin
		ap[i,*]=ap[first,*]+ap1[i]-ap[first,np/2]
	endfor
endif

if last lt nap1-1 then begin
	for i=last+1,nap1-1 do begin
		ap[i,*]=ap[last,*]+ap1[i]-ap[last,np/2]
	endfor
endif


end
