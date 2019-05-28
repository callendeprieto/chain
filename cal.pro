pro cal,i,f1,f2

;load,'thar.ascii',t
;load,'table.dat',l
;t=double(t)
;l=double(l)
t=readfits('thar.fits')
l=readfits('table.fits')
rs,'ce.fits',ce,cev,w=x

loadct,12
;for i=1,n_elements(x[*,0])-1 do begin
  loci=-1
  y=-1
  loci2=-1
  y2=-1
  plot,x[i,*],ce[i,*],charsi=2
  peaks,ce[i,*],loci,y,fwhm=4.,thres=max(ce[i,*])/f1
  w=where(t[0,*] gt min(x[i,*]) and t[0,*] lt max(x[i,*]))
  scale=max(t[1,w])/max(ce[i,*])
  tt0=t[0,w]
  tt=smooth(t[1,w]/scale,30)
  oplot,tt0,tt,col=180
  peaks,tt,loci2,y2,fwhm=4.0,thres=max(ce[i,*])/f2
  loci=loci[1:n_elements(loci)-1]
  loci2=loci2[1:n_elements(loci2)-1]
  print,n_elements(loci),'peaks in hors'
  print,'at locations',loci
  for j=0,n_elements(loci)-1 do begin
      w1=x[i,loci[j]]
      plots,[w1,w1],[max(ce[i,*])*1.05,max(ce[i,*])*1.15]
  endfor
  print,n_elements(loci2),'peaks in mcd'
  print,'at wavelengths ... '
  for j=0,n_elements(loci2)-1 do begin
    w0=tt0[loci2[j]]
    ww=where(abs(l[0,*]-w0) eq min(abs(l[0,*]-w0)))
    w1=l[0,ww]
    plots,[w1,w1],[max(ce[i,*])*1.1,max(ce[i,*])*1.2],col=180
    print,w0,w1
  endfor
  ;stop
;endfor

end

