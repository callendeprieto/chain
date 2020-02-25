pro xm,xfile

;
; Continuum normalization and order merging
; It uses an extracted flat present in the working dir (xflat.fits)
;
; IN: xfile -- string Name of an x*fits file (extracted chain spectrum)
;

rs,xfile,y,v,w=x,hd=hd
norder=n_elements(y[*,0])

rs,'xflat.fits',yf,vf,w=xf
margin=50
yy=y[*,margin:2048-margin]
yyf=yf[*,margin:2048-margin]
vv=v[*,margin:2048-margin]
xx=x[*,margin:2048-margin]
for i=0,norder-1 do yyf[i,*]=yyf[i,*]/mean(yyf[i,*])
yy2=yy/yyf
vv2=vv/yyf^2
for i=0,norder-1 do begin
  plot,yy2[i,*]
  ;c=cmedian(yy2[i,*],300,percentile=0.85)
  ;coef=poly_fit(findgen(n_elements(c)),c,8,yfit=c2)
  ;oplot,c,col=100
  ;oplot,c2,col=180
  if i lt 14 then corder=6 else corder=9
  continuum,corder,10,1.,5.,yy2[i,*],c3
  oplot,c3,col=180,thick=3
  ;stop
  yy2[i,*]=yy2[i,*]/c3
  vv2[i,*]=vv2[i,*]/c3^2
endfor

cmerge,xx,yy2,vv2,x2,y2,v2
x2=x2[margin:n_elements(x2)-1]
y2=y2[margin:n_elements(y2)-1]
v2=v2[margin:n_elements(v2)-1]
wout=where(y2 lt 0.0)
if max(wout) gt -1 then begin
  y2[wout]=0.0
endif

data=replicate({wavelength:x2[0],flux:y2[0],var:v2[0]},n_elements(x2))
data.wavelength=x2
data.flux=y2
data.var=v2

mwrfits,data,'m'+strmid(xfile,1),hd[7:n_elements(hd)-1],/create

end
