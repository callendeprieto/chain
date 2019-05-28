pro	map,s,t,fwhm=fwhm,w=w

;
;	Identifies features in an extracted spectrum
;
;	IN:    s    - float array		Spectrum
;			np x nap
;
;	OUT: 	t   - float array		Table with identifed features
;			4 x nfeatures		aperture, pixel, flux, fwhm
;
;	KEYWORDS: fwhm - float	expected fwhm of features (peaks); default is 3
;		  w	- float array		Wavelegnths. If provided, then 
;						the return values in table t for pixel
;						and fwhm will be in the same units as w
;


if n_elements(t) gt 0 then undefine,t
if n_elements(fwhm) eq 0 then fwhm=3.0
sigma_to_fwhm=2.0*sqrt(-2.0*alog(0.5))

ss=size(s)
if ss[0] ne 2 then begin
	print,'% MAP: error the input array is not a 2D spectrum'
	return
endif

np=ss[2]
nap=ss[1]

for i=0,nap-1 do begin

        plot,s[i,*]
	peaks,s[i,*],loc,val,thres=median(s[i,*])*2.+1000.,/pl
        wait,5

	nloc=n_elements(loc)
	print,' aperture=',i,' --- ',n_elements(loc),' features'

	loc=array(loc)
	val=array(val)
	if nloc gt 0 then begin
		ii=replicate(i,nloc)
		width=fltarr(nloc)
		for j=0,nloc-1 do begin
			jrange=indgen(5*fwhm)-fix(2.5*fwhm)+loc[j]
			g=gaussfit(jrange,s[i,jrange],nter=4,aa)
			width[j]=aa[2]*sigma_to_fwhm
			if n_elements(w) gt 0 then begin 
			; convert width to units of wavelength
			  ds=w[i,jrange]-shift(w[i,jrange],1)
			  ds=ds[1:n_elements(ss)-1]
			  width[j]=width[j]*mean(ds)
			endif
			;plot,jrange,s[jrange,i],psym=10
			;oplot,jrange,g,col=180
			;print,aa
			;stop
		endfor
	endif

	if n_elements(w) gt 0 then loc=interpol(w[i,i],dindgen(np),loc,/spl)

	if n_elements(t) gt 0 and n_elements(loc) gt 0 then $
                t=[[t],[transpose([[ii],[loc],[val],[width]])]]

	if n_elements(t) eq 0 and n_elements(loc) gt 0 then $
		t=transpose([[ii],[loc],[val],[width]])

endfor

end

