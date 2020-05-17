pro calrefine,cc,ss

;+
; Fit polynomials to the wavelength calibration coefficients to 
; enforce smooth variation in the cross-dispersion direction
;
; IN: cc -- array with the wavelength calibration polynomial coefficients (norder,ncoef)
;     ss -- array with the uncertainties in the calibration coefficients
;     
; OUT: cc and ss are updated on output
;-


;fit the coefficients as a function of order
ncoef=n_elements(cc[0,*])
norder=n_elements(cc[*,0])
;!p.multi=[0,3,ceil(ncoef/3.)
degree=replicate(0,ncoef)
degree[0]=8
degree[1]=2
!p.multi=[0,2,ncoef]
enhancement_errors=1.
for i=0,ncoef-1 do begin
  ploterror,cc[*,i],ss[*,i]*enhancement_errors,psy=-4,charsi=2
  xcross=dindgen(norder)
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
endfor



end
