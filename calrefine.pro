pro calrefine,cc,ss,plot=plot

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
degree[0]=4
degree[1]=2
enhancement_errors=1.
if keyword_set(plot) then begin
  loadct,12
  !p.multi=[0,2,ncoef]
endif
for i=0,ncoef-1 do begin

  xcross=dindgen(norder)
  
  if keyword_set(plot) then begin

    plot,xcross,cc[*,i],psym=-4,xtit='Aperture #',$
      title='Coefficient #'+string(i),charsi=1.6
    oploterror,xcross,cc[*,i],replicate(0.0,norder),$
     ss[*,i]*enhancement_errors,psym=4,charsi=2

  endif 

  coef=poly_fit(xcross,cc[*,i],degree[i],yfit=yfit,yerror=yerror,sigma=sigma,$
	measure_errors=ss[*,i]*enhancement_errors)
  coef2=poly_fit(xcross,cc[*,i],degree[i]+1,yfit=yfit2,yerror=yerror2,sigma=sigma2,$
	measure_errors=ss[*,i]*enhancement_errors)

  if n_elements(yfit) eq 1 then yfit=replicate(yfit,norder)

  if keyword_set(plot) then begin
    oplot,xcross,yfit,col=180,thick=3
    plot,xcross,yfit-cc[*,i],psy=-4,yr=[-yerror,yerror]/3.,charsi=1.6,$
      xtit='Aperture #',ytit='Residuals'
    oplot,xcross,yfit2-cc[*,i],psy=-2,col=100
    xyouts,max(xcross)*0.7,yerror/3.*0.7,'order '+string(degree[i])+ '(adopted)',charsi=2
    xyouts,max(xcross)*0.7,yerror/3.*0.3,'order '+string(degree[i]+1),col=100,charsi=2
  endif

  print,'mean-std-error(order='+strcompress(string(degree[i]),/rem)+')-error(order='+strcompress(string(degree[i]+1),/rem)+')',mean(cc[*,i]),stddev(cc[*,i]),yerror,yerror2
  w=where(ss[*,i] lt median(ss[*,i]))
  ;coco=ladfit(xcross[w],cc[w,i])
  ;oplot,xcross,poly(xcross,coco)-cc[*,i],col=60
  print,'aperture=',i,' deg=',degree[i],' mad(residuals)=',mad(yfit-cc[*,i])
  print,'aperture=',i,' deg=',degree[i]+1,' mad(residuals)=',mad(yfit2-cc[*,i])
  ;stop
  ss[*,i]=replicate(mad(yfit-cc[*,i]),norder)
  cc[*,i]=poly(xcross,coef)
endfor



end
