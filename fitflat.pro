pro  fitflat,flat,fit


order=6
rs,flat,f

sf=size(f)
if sf[0] ne 2 then begin
	message,'Dimension of data array must be 2!'
	return
endif

npix=sf[2]
norder=sf[1]
fit=dblarr(norder,npix)

x=findgen(npix)
for i=0,norder-1 do begin
  c=poly_fit(x,cmedian(f[i,*],100),order)
  plot,x,f[i,*]
  fit[i,*]=poly(x,c)
  oplot,x,fit[i,*],col=180,thick=2
endfor

end
