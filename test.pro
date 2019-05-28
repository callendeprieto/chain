b=readfits('0001610105-20180629-HORS-Bias.fits',bhd)         ; bias
s=readfits('0001610254-20180629-HORS-Spectroscopy.fits',shd) ; skylight
c=readfits('0001610255-20180629-HORS-Cal.fits',chd)          ; thar
h=readfits('0001610354-20180629-HORS-Spectroscopy.fits',ohd) ; hd140283

dispdir,s,idisp  ;find out the dispersion direction from a sky exposure
s=s-b            ;subtract bias frame
c=c-b
h=h-b
cut,s            ;remove overscan part of the detector
cut,c
cut,h

gain=sxpar(shd,'GAIN')
rdnoise=sxpar(shd,'RDNOISE')

sv=s*gain+rdnoise^2  ;compute variance images
cv=c*gain+rdnoise^2
hv=h*gain+rdnoise^2

find1,s,idisp,ap1,delta1  ;find positions and width of orders
                          ;assuming disp. direction is aligned with CCD
;collapse1,s,idisp,ap1,ys
;collapse1,c,idisp,ap1,yc
;collapse1,h,idisp,ap1,yh

find,s,idisp,ap1,ap       ;determine curvature of orders

inspect,s,idisp,ap,delta1

collapse,s,idisp,ap,delta1,se,vf=sv,vs=sev ; extract
collapse,c,idisp,ap,delta1,ce,vf=cv,vs=cev
collapse,h,idisp,ap,delta1,he,vf=hv,vs=hev

;rlambda,n_elements(se[0,*]),x  ; very rough wavelength cal.

xcal,ce,x ; wavelength calibration


;window,0
;inspect1,s[3000:4000,*],idisp,ap1-3000.
;window,2
;inspect,s[3000:4000,*],idisp,ap-3000.,delta1


end

