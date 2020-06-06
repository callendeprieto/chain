pro mkthar

;builds an horus thar file from an approximately calibrated HORuS Thar
;which is compared with the McDonald ThAr atlas

threshold=0.01 ; relative height of peaks to consider
epsilon=0.0005 ; pixel space error for matching features

openw,11,'horus_thar1.dat'
rs,'x0002364843-20191120-HORS-Cal.fits',yy,w=xx
window,0,xsize=900,ysize=600
loadct,12
norder=27
for i=0,norder-1 do begin
  x=transpose(xx[i,*])*10.
  y=transpose(yy[i,*])*10.
  yc=median(y,100)
  !p.multi=[0,2,2]
  read_thar,min(x),max(x),x2,y2;,/pl
  gconv,x2,y2,3,nx2,ny2,/ori
  y2=ny2
  x2=nx2  
  y2=interpol(y2,x2,x,/spli)
  ;oplot,x,y,col=180
  ;oplot,x,yc,col=100
  y2=y2/max(y2)
  y=y-yc
  y=y/max(y)

  print,'aperture=',i

  ;plot,x,y2
  ;oplot,x,y,col=180


  ;first determine the appropriate fwhm
  peaks,y2,l2,h2,fwhm=3.0,threshold=threshold
  ws=l2[sort(h2)]
  wi=floor(ws[n_elements(h2)-3])
  w=indgen(20)+wi-10
  g=gaussfit(x[w],y2[w],nter=3,a)
  fwhm=2.355*a[2]/median(x-shift(x,1))

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

  ;now crossmatch peaks
  match,r,r2,wl,wl2,epsil=epsilon

  plot,y2,yr=[0.0,1.3],/ystyle,xtit='pixel',ytit='flux'
  oplot,y,col=180
  for j=0,n_elements(wl2)-1 do oplot,[l2[wl2[j]],l2[wl2[j]]],[1.02,1.07]
  for j=0,n_elements(wl)-1 do oplot,[l[wl[j]],l[wl[j]]],[1.1,1.15],col=180
  xyouts,n_elements(y2)*0.8,max(y2)*0.7,'McDonald',charsi=2
  xyouts,n_elements(y2)*0.8,max(y2)*0.8,'HORuS',col=180,charsi=2
  ;help,wl2,wl
  print,n_elements(wl2),' peaks cross-matched'

  l2=l2[wl2]
  l=l[wl]

  ref=interpol(x,dindgen(n_elements(x)),l2)
  plot,l,ref,psy=2,yr=[min(ref),max(ref)],xtit='line center (measured pixel)',ytit='assigned line center wavelength (from McD atlas)'
  coef=poly_fit(l,ref,2,yfit=yfit,yerror=yerror,sigma=sigma)
  oplot,l,poly(l,coef),col=180
  plot,l,(poly(l,coef)-ref)/ref,yr=[-yerror/mean(ref)*3,yerror/mean(ref)*3],charsi=1.6,psy=4,xtit='line center pixel',ytitle='wavelength residual',title='2nd order fit'
  mmad=mad(poly(l,coef)-ref)
  print,'order,rms,mad,disp=',i,yerror,mmad,coef[1]

  ;cleanup and repeat
  w=where(poly(l,coef)-ref lt mmad)
  l=l[w]
  l2=l2[w]

  ref=interpol(x,dindgen(n_elements(x)),l2)
  coef=poly_fit(l,ref,4,yfit=yfit,yerror=yerror,sigma=sigma)
  plot,l,(poly(l,coef)-ref)/ref,yr=[-yerror/mean(ref)*3,yerror/mean(ref)*3],charsi=1.6,psy=4,xtit='line center pixel',ytitle='wavelength residual',title='4th order fit'
  mmad=mad(poly(l,coef)-ref)
  print,'order,rms,mad,disp=',i,yerror,mmad,coef[1]
  printf,11,i+1,n_elements(l),4
  printf,11,transpose(coef)
  printf,11,sigma
  for j=0,n_elements(l)-1 do printf,11,i+1,j+1,l[j],ref[j]/10.d0

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
close,11
ss=transpose(ss)

tharrefine,'horus_thar1.dat','horus_thar.dat'

end
