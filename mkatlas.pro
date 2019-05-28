rs,'ce.fits',ce,cev,w=x

openr,11,'hors_thar.dat'
;psinit,/vector,/landscape
for i=0,26 do begin

  plot,x[i,*],ce[i,*],chars=1.8,xtitle='Wavelength (nm)',$
	ytitle='counts',title='ThAr',yr=[0.0,max(ce[i,*])*1.3]

  i2=0
  nlines=0
  order=0
  readf,11,i2,nlines,order
  coef=dblarr(order+1)
  readf,11,coef
  xx=poly(dindgen(4088),coef)

  k=0
  for j=0,nlines-1 do begin

    i2=0
    j2=0
    pix=0.0d0
    l0=0.0d0
    readf,11,i2,j2,pix,l0
    wline=where(x[i,*] ge l0-0.05 and x[i,*] le l0+0.05)
    flux=mean(ce[i,wline])
    if flux ge max(ce[i,*])*0.05 then begin
	    oplot,[l0,l0],[max(ce[i,*])*1.1,max(ce[i,*])*1.2]
            print,l0*10.,format='(f10.5)"
	    xyouts,l0*0.9975,max(ce[i,*])*1.1*(1.-k/10.),string(l0),charsi=1.3
            print,1.-k/10.
            if (1.-k/10.) gt 0.3 then k=k+1 else k=k-8
    endif

  endfor
  write_png,'horus_thar'+string(i,format='(i02)')+'.png',tvrd(/true)
endfor
close,11
;psterm,/noplot,file='horus_thar.ps'
end