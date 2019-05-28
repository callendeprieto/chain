openr,1,'caldata.dat'
openw,2,'caldata2.dat'
loadct,12
while not eof(1) do begin
  readf,1,ap
  readf,1,ndata
  order=1
  pix=dblarr(ndata)
  lam=dblarr(ndata)
  readf,1,pix
  readf,1,lam
  plot,pix,lam,psy=-2,charsi=2,xr=[min(pix)-10,max(pix)+10],$
	yr=[min(lam)-10.,max(lam)+10.],/ystyl,/xstyl,title=string(ap)
  c=poly_fit(pix,lam,order)

  print,'ap#/ndata/order/c/mean/std=',ap,ndata,order,c,mean(poly(pix,c)-lam),stddev(poly(pix,c)-lam)
  oplot,pix,poly(pix,c),col=180,thick=2,psy=-4,symsize=2
  ;plot,pix,poly(pix,c)-lam,psy=4
  printf,2,ap,c,stddev(poly(pix,c)-lam),format='(i2,f10.4,f10.4,f10.4)'
endwhile
close,1
close,2

end

