pro tharrefine,tharfile,tharfile2

;+
; Fit polynomials to the wavelength calibration coefficients to 
; enforce smooth variation in the cross-dispersion direction
;
; IN: tharfile -- file with wavelength calibration polynomials 
;     tarfile2 -- file for the updated data
;
;-

openr,11,tharfile
while not eof(11) do begin
 i=0
 nlines=0
 order=0
 readf,11,i,nlines,order
 coef=dblarr(order+1)
 sigma=dblarr(order+1)
 pix=dblarr(nlines)
 lambda=dblarr(nlines)
 readf,11,coef
 readf,11,sigma
 for j=0,nlines-1 do readf,11,i2,jj,pix[j],lambda[j]
 if i eq 1 then begin
   order0=order
   cc=coef
   ss=sigma
 endif else begin
  if order ne order0 then stop,'the polynomial degree must be the same for all the orders'
  cc=[[cc],[coef]]
  ss=[[ss],[sigma]]
 endelse
endwhile
close,11

cc=transpose(cc)
ss=transpose(ss)

;fit the coefficients as a function of order
ncoef=n_elements(cc[0,*])
norder=n_elements(cc[*,0])
;!p.multi=[0,3,ceil(ncoef/3.)
degree=replicate(0,ncoef)
degree[0]=8
degree[1]=2
!p.multi=[0,1,2]
enhancement_errors=1.
for i=0,ncoef-1 do begin
  ploterror,cc[*,i],ss[*,i]*enhancement_errors,psy=-4,charsi=2
  xcross=dindgen(norder)
  stop
  coef=poly_fit(xcross,cc[*,i],degree[i],yfit=yfit,yerror=yerror,sigma=sigma,$
	measure_errors=ss[*,i]*enhancement_errors)
  coef2=poly_fit(xcross,cc[*,i],degree[i]+1,yfit=yfit2,yerror=yerror2,sigma=sigma2,$
	measure_errors=ss[*,i]*enhancement_errors)

  if n_elements(yfit) eq 1 then yfit=replicate(yfit,ncoef)
  oplot,xcross,yfit,col=180,thick=3
  plot,xcross,yfit-cc[*,i],psy=-4,yr=[-yerror,yerror]/3.,charsi=2
  oplot,xcross,yfit2-cc[*,i],psy=-2,col=100
  print,mean(cc[*,i]),stddev(cc[*,i]),yerror,yerror2
  w=where(ss[*,i] lt median(ss[*,i]))
  ;coco=ladfit(xcross[w],cc[w,i])
  ;oplot,xcross,poly(xcross,coco)-cc[*,i],col=60
  print,'aperture=',i,' deg=',degree[i],' error=',mad(yfit-cc[*,i])
  print,'aperture=',i,' deg=',degree[i]+1,' error=',mad(yfit2-cc[*,i])
  ;stop
  ss[*,i]=replicate(mad(yfit-cc[*,i]),norder)
  cc[*,i]=poly(xcross,coef)
  stop
endfor

;write tharfile2
openr,11,tharfile
openw,12,tharfile2
while not eof(11) do begin
 i=0
 nlines=0
 order=0
 readf,11,i,nlines,order
 printf,12,i,nlines,order
 coef=dblarr(order+1)
 sigma=dblarr(order+1)
 readf,11,coef
 printf,12,transpose(cc[i-1,*])
 readf,11,sigma
 printf,12,transpose(ss[i-1,*])
 for j=0,nlines-1 do begin
   i2=0
   jj=0
   pix=0.d0
   lambda=0.0d0
   readf,11,i2,jj,pix,lambda
   printf,12,i2,jj,pix,poly(pix,cc[i-1,*])/10.d0
 endfor
endwhile
close,11
close,12


end
