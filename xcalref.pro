rs,'ce3.fits',ce,cev,w=x

xnew=x

t=mrdfits('thar.fits',0,hd)
xt=transpose(t[0,*])/10.
yt=transpose(t[1,*])
gconv,xt,yt,6.,xt2,yt2
xt=xt2
yt=yt2

load,'table.dat',tab,skip=1
tab=double(tab)
tab[0,*]=tab[0,*]/10.

loadct,12

openw,10,'results3.dat'
openw,11,'hors_thar.dat'
for i=0,26 do begin

  ptol=3.2 ; pixels from the expected location to accept a line fit
  tol=0.07 ; nm around nominal central wavelength to search for a line in the HORS spectrum

  plot,x[i,*],ce[i,*]
  oplot,xt,yt*5.+median(ce[i,*]),col=180

  wlines=where(tab[0,*] ge min(x[i,*])+.5 and tab[0,*] le max(x[i,*])-.5 and tab[1,*] lt 0.001)

  print,'aperture = ',i+1
  print,'found '+string(n_elements(wlines))+' candidate calibration lines'
  for j=0,n_elements(wlines)-1 do $
    oplot,replicate(tab[0,wlines[j]],2),[max(ce[i,*])/2.,max(ce[i,*])]

  ;lets measure
  ngood=0
  pixels=-1.d0
  lambdas=-1.d0
  for j=0,n_elements(wlines)-1 do begin
   
     ;first McD atlas
     l0=tab[0,wlines[j]]
     wxt=where(abs(xt-l0) lt 0.05)
     wx=where(abs(x[i,*]-l0) lt tol)

     g1=gaussfit(dindgen(n_elements(wxt)),yt[wxt],nter=4,a1)
     !p.multi=[0,1,2]
     plot,yt[wxt],title=string(a1[0])+string(a1[1])+string(a1[2])
     oplot,g1,col=180

     ;select only those that fit well in the atlas
     nwxt=n_elements(wxt)
     if a1[0] gt 100. and a1[1] gt nwxt/2.-1.-1.5 and a1[1] lt nwxt/2.-1.+1.5 and abs(a1[2]) gt 1.0 and abs(a1[2]) lt 1.6 then begin

       ;now fit the HORS data
       g=gaussfit(dindgen(n_elements(wx)),ce[i,wx],nter=4,a)

       nwx=n_elements(wx)
       if a[0] gt 100. and a[1] gt nwx/2.-1.-ptol and a[1] lt nwx/2.-1.+ptol and abs(a[2]) gt 1.0 and abs(a[2]) lt 3.3 then begin


       plot,findgen(nwx)*nwxt/nwx,ce[i,wx],$
	 title=string(a[0])+string(a[1])+string(a[2]),psym=-2
       oplot,findgen(nwx)*nwxt/nwx,g,col=180
       ngood=ngood+1
       pixels=[pixels,a[1]+wx[0]]
       lambdas=[lambdas,l0]

       ;wait,2

       endif
     endif

  endfor
 
  print,ngood,' successfully measured lines'
  pixels=pixels[1:ngood-1]
  lambdas=lambdas[1:ngood-1]

  c1=poly_fit(pixels,lambdas,1,yerror=yerror1)
  c2=poly_fit(pixels,lambdas,2,yerror=yerror2)
  c3=poly_fit(pixels,lambdas,3,yerror=yerror3)  
  c4=poly_fit(pixels,lambdas,4,yerror=yerror4)
  print,i+1,n_elements(pixels),median(x[i,*]-shift(x[i,*],1)),$
	yerror1,yerror2,yerror3,yerror4

  wgood=where(abs(lambdas-poly(pixels,c4)) lt 0.003)
  pixels=pixels[wgood]
  lambdas=lambdas[wgood]
  c1=poly_fit(pixels,lambdas,1,yerror=yerror1)
  c2=poly_fit(pixels,lambdas,2,yerror=yerror2)
  c3=poly_fit(pixels,lambdas,3,yerror=yerror3)
  c4=poly_fit(pixels,lambdas,4,yerror=yerror4)
  print,i+1,n_elements(pixels),median(x[i,*]-shift(x[i,*],1)),$
	yerror1,yerror2,yerror3,yerror4

  xnew[i,*]=poly(dindgen(4088),c4)

  printf,10,i+1,n_elements(pixels),median(xnew[i,*]-shift(xnew[i,*],1)),$
	yerror1,yerror2,yerror3,yerror4,format='(i3,x,i3,5(1x,f14.7))'

  printf,11,i+1,n_elements(pixels),4,format='(i3,x,i3,x,i3)'
  printf,11,c4,format='(5(x,e20.9))'
  for j=0,n_elements(pixels)-1 do begin
    printf,11,i+1,j+1,pixels[j],lambdas[j],format='(i3,x,i3,2(x,f20.9))'
  endfor

  !p.multi=[0,0,0]
  plot,pixels,lambdas,psy=-2
  plot,pixels,lambdas-poly(pixels,c4),yr=[-0.005,0.005],psy=2
  oplot,[0,1e4],[-0.002,-0.002],linestyle=2
  oplot,[0,1e4],[0.002,0.002],linestyle=2
  ;wait,3



endfor
close,10
close,11
end

