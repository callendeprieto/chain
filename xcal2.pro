pro xcal2,c,x,calstats=calstats,chaindir=chaindir,plot=plot,cc=cc,ss=ss

;
;	IN:  c -float array 	Collapsed ThAr spectrum
;	OUT: x -double array 	Wavelengths corresponding to flux array c
;
;	KEYWORDS: calstats   7x27 array with 7 columns giving: aperture, nlines,
 
;		               dispersion, yerror1, yerror2, yerror3, yerror4
;
;		   where yerrorN is the uncertainty for a N-th order polynomial calib.
;
;		chaindir - string  Path to the chain (both software and reference files)
;
;		cc - fltarr  Array with polynomial coefficients  (26 apertures,5 coeffs)
;		ss - fltarr  Array with uncertainties in the coefficients cc (26x5)
;


threshold=0.0001 ; relative height of peaks to consider
epsilon=0.0001 ; pixel space error for matching features
min_nlines=50  ; minimum number of lines to cross-match 
               ; epsilon will be increased until reaching it)

home=getenv('HOME')
if not keyword_set(chaindir) then chaindir=home+'/idl/chain'


norder=27
nx=n_elements(c[0,*])
x=dblarr(norder,nx)
calstats=dblarr(7,norder)

;read reference spectrum (used only to find wavelength range for each order)
rs,chaindir+'/x0002364843-20191120-HORS-Cal.fits',ydummy,w=xx 
;coe=poly_fit(findgen(27),xx[*,2047]-xx[*,0],2)
;range=poly(findgen(27),coe)


if keyword_set(plot) then begin
  window,0,xsize=900,ysize=600
  loadct,12
endif

for i=0,norder-1 do begin
  x1=transpose(xx[i,*])*10.
  ;x1=range[i]*10./2048.*(dindgen(2048)-1024.)+mean(xx[i,*])*10.
  y1=transpose(c[i,*])
  yc=median(y1,100)
  !p.multi=[0,2,2]
  read_thar,min(x1),max(x1),x2,y2;,/pl
  ;read_thar,mean(x1)-range[i]*10./2.,mean(x1)+range[i]*10./2.,x2,y2;,/pl
  gconv,x2,y2,3,nx2,ny2,/ori
  y2=ny2
  x2=nx2  
  y2=interpol(y2,x2,x1,/spli)
  ;oplot,x,y,col=180
  ;oplot,x,yc,col=100
  y2=y2/max(y2)
  y=y1-yc
  y=y/max(y)

  print,'aperture=',i


  ;first determine the appropriate fwhm
  peaks,y2,l2,h2,fwhm=3.0,threshold=threshold
  ws=l2[sort(h2)]
  wi=floor(ws[n_elements(h2)-3])
  w=indgen(20)+wi-10
  g=gaussfit(x1[w],y2[w],nter=3,a)
  fwhm=2.355*a[2]/median(x1-shift(x1,1))

  ;now measure peaks
  peaks,y2,l2,h2,fwhm=fwhm,threshold=threshold
  peaks,y,l,h,fwhm=fwhm,threshold=threshold
  ;help,l,l2
  print,n_elements(l),n_elements(l2),' peaks found in the HORuS / McD spectra'
  ;to coords rel. to highest peak
  ws=l2[sort(h2)]
  lm=n_elements(y)/2.
  r2=(l2-lm)/lm
  ws=l[sort(h)]
  r=(l-lm)/lm


  epsil=epsilon
  wl2=0

  while n_elements(wl2) lt min_nlines do begin 
    ;now crossmatch peaks
    match,r,r2,wl,wl2,epsil=epsil
    print,n_elements(wl2),' peaks cross-matched'
    epsil=epsil*2.
  endwhile

  if keyword_set(plot) then begin
    plot,y2,yr=[0.0,1.3],/ystyle,xtit='pixel',ytit='flux',title='aperture #'+string(i)
    oplot,y,col=180
    for j=0,n_elements(wl2)-1 do oplot,[l2[wl2[j]],l2[wl2[j]]],[1.02,1.07]
    for j=0,n_elements(wl)-1 do oplot,[l[wl[j]],l[wl[j]]],[1.1,1.15],col=180
    xyouts,n_elements(y2)*0.8,max(y2)*0.7,'McDonald',charsi=2
    xyouts,n_elements(y2)*0.8,max(y2)*0.8,'HORuS',col=180,charsi=2
    ;help,wl2,wl
  endif

  l2=l2[wl2]
  l=l[wl]

  ref=interpol(x1,dindgen(n_elements(x1)),l2)
  coef=poly_fit(l,ref,2,yfit=yfit,yerror=yerror,sigma=sigma)

  if keyword_set(plot) then begin
    plot,l,ref,psy=2,yr=[min(ref),max(ref)],xtit='line center (measured pixel)',ytit='assigned line center wavelength (from McD atlas)'
    oplot,l,poly(l,coef),col=180
    plot,l,(poly(l,coef)-ref)/ref,yr=[-yerror/mean(ref)*3,yerror/mean(ref)*3],charsi=1.6,psy=4,xtit='line center pixel',ytitle='wavelength residual',title='2nd order fit'
  endif

  mmad=mad(poly(l,coef)-ref)
  print,'order,rms,mad,disp=',i,yerror,mmad,coef[1]

  ;cleanup and repeat
  w=where(poly(l,coef)-ref lt mmad)
  l=l[w]
  l2=l2[w]

  ref=interpol(x1,dindgen(n_elements(x1)),l2)
  coef1=poly_fit(l,ref,1,yfit=yfit1,yerror=yerror1,sigma=sigma1)
  coef2=poly_fit(l,ref,2,yfit=yfit2,yerror=yerror2,sigma=sigma2)
  coef3=poly_fit(l,ref,3,yfit=yfit3,yerror=yerror3,sigma=sigma3)
  coef=poly_fit(l,ref,4,yfit=yfit,yerror=yerror,sigma=sigma)

  if keyword_set(plot) then begin
    plot,l,(poly(l,coef)-ref)/ref,yr=[-yerror/mean(ref)*3,yerror/mean(ref)*3],charsi=1.6,psy=4,xtit='line center pixel',ytitle='wavelength residual',title='4th order fit'
  endif

  mmad=mad(poly(l,coef)-ref)
  print,'order,rms,mad,disp=',i,yerror,mmad,coef[1]

  x[i,*]=poly(dindgen(nx),coef)/10.

  ;fill-in calstats
  calstats[*,i]=[i,n_elements(w),median(x[i,*]-shift(x[i,*],1)),$
	yerror1/10.,yerror2/10.,yerror3/10.,yerror/10.]


  ;store coefficients in array
  if i eq 0 then begin
    cc=coef 
    ss=sigma
  endif else begin
    cc=[cc,[coef]]
    ss=[[ss],[sigma]]
  endelse
  ;stop

endfor
ss=transpose(ss)

;enforce coherence in the calibration across orders
;calrefine,cc,ss,/plot
;
;for i=0,26 do begin
;  x[i,*]=poly(dindgen(nx),cc[i,*])
;endfor


end
