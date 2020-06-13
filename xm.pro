pro xm,xfile,xmformat=xmformat

;
; Continuum normalization and order merging
; It uses an extracted flat present in the working dir (xflat.fits)
;
; IN: xfile -- string Name of an x*fits file (extracted chain spectrum)
;

if not keyword_set(xmformat) then xmformat='vo'

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
  ;plot,yy2[i,*]
  ;;c=cmedian(yy2[i,*],300,percentile=0.85)
  ;;coef=poly_fit(findgen(n_elements(c)),c,8,yfit=c2)
  ;;oplot,c,col=100
  ;;oplot,c2,col=180
  if i lt 14 then corder=6 else corder=9
  continuum,corder,10,1.,5.,yy2[i,*],c3
  ;oplot,c3,col=180,thick=3
  ;;;stop
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


case xmformat of
  'vo': begin
             data=replicate({wavelength:x2[0],flux:y2[0],var:v2[0]},n_elements(x2))
             data.wavelength=x2
             data.flux=y2
             data.var=v2
             mwrfits,data,'v'+strmid(xfile,1),hd[7:n_elements(hd)-1],/create
        end
   'three-col': begin
             data2=dblarr(3,n_elements(x2))
             data2[0,*]=x2
             data2[1,*]=y2
             data2[2,*]=v2
             mwrfits,data2,'t'+strmid(xfile,1),hd[7:n_elements(hd)-1],/create
                end
    'rana': begin
             delta=5d-3
             nfreq=(max(x2)-min(x2))/delta
             xx=min(x2)+dindgen(nfreq)*delta
             yy=interpol(y2,x2,xx)
             sxaddpar,hd,'crpix1',1
             sxaddpar,hd,'crval1',min(xx)*10.d0 ; AA
             sxaddpar,hd,'cdelt1',delta*10.d0   ; AA
             writefits,'r'+strmid(xfile,1),yy,hd
            end
     else: begin
             print,'% XM: xmformat not recognized -- valid options are:'
             print,'% XM:      vo, three-col, rana'
           end
endcase


end
