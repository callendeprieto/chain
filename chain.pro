pro chain,bin=bin,logfile=logfile

if not keyword_set(logfile) then logfile='logfile'
openw,10,logfile

nboost=3 ;factor to increase spatial sampling when binning is on

;make data inventory 
inventory,st,bin=bin

print,'make inventory of fits files ...'
print,'                           filename                object  obstype    mjd0     exptime [s]'
printf,10,'make inventory of fits files ...'
printf,10,'                           filename                object  obstype    mjd0     exptime [s]'
for i=0,n_elements(st.obstype)-1 do begin
	printf,10,st[i].filename,st[i].object,strmid(st[i].obstype,0,3),$
		st[i].mjd0,st[i].exptime,$
		format='(x,a45,x,a12,x,a3,x,f12.5,x,f7.0)'
	print,st[i].filename,st[i].object,strmid(st[i].obstype,0,3),$
		st[i].mjd0,st[i].exptime,$
		format='(x,a45,x,a12,x,a3,x,f12.5,x,f7.0)'
endfor


wcal=where(strmid(st.obstype,0,3) eq 'Cal')
wspe=where(strmid(st.obstype,0,3) eq 'Spe')
wfla=where(strmid(st.obstype,0,3) eq 'Fla')
wob=[wcal,wspe]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;if keyword_set(bin) then $
;	;s=readfits('0000000000-20190312-HORS-Spectroscopy.fits',shd) else $
;	s=readfits('0002304935-20190926-HORS-Spectroscopy.fits',shd) else $
;	s=readfits('0001947116-20190312-HORS-Spectroscopy.fits',shd) ; skylight


;find order information from skylight spectrum
;print,'find order information from skylight spectrum...'
;printf,10,'find order information from skylight spectrum...'
;dispdir,s,idisp  ;find out the dispersion direction from a sky exposure


;find positions and width of orders
;assuming disp. direction is aligned with CCD
;hbias,s
;if keyword_set(bin) then begin
;        s=rebin(s,514*nboost,2048) ; icreasing the sampling in the spatial dir. 
;	hfind, s, idisp, ap1, delta1, fwhm=4*nboost, smoothinglength=4*nboost, edge=19 ; 19.3
;endif else begin
;	hfind, s, idisp, ap1, delta1, fwhm=10.0, smoothinglength=55.
;endelse


;if n_elements(ap1) ne 27 then begin
;	print,'% CHAIN: ERROR -- did not find exactly 27 orders!'
;	print,'               -- this is a requirement for xcal.'
;	printf,10,'% CHAIN: ERROR -- did not find exactly 27 orders!'
;	printf,10,'               -- this is a requirement for xcal.'
;	close,10
;       stop
;endif


;determine curvature of orders
;traze, s, idisp, ap1, ap , smoothinglength=8., order=1,  width=0.7, last=21
;inspect, s, idisp, ap, delta1
;window,2
;inspect, s[0:200*nboost,*], idisp, ap, delta1
;window,3
;inspect, s[200*nboost:300*nboost,*], idisp, ap-200*nboost, delta1
;window,4
;help,s
;print,fix(n_elements(s[*,0])*0.9),n_elements(s[0,*])
;inspect, s[300*nboost:fix(n_elements(s[*,0])*0.9),*], idisp, ap-300*nboost, delta1

ap1=readfits('/home/callende/idl/hors/rap1.fits')
ap=readfits('/home/callende/idl/hors/rap.fits')
delta1=13.


;average flats and extract average
printf,10,'creating average flat ...'
print,'creating average flat ...'
imadd,st[wfla].filename,'flat.fits',/av
f = readfits('flat.fits',header)


;derive spatial-direction offset of orders using flat
rflat=readfits('/home/callende/idl/hors/rflat.fits')
trflat=total(rflat,2)
tflat=total(f,2)
if not keyword_set(bin) then tflat=rebin(tflat,514)
xc,trflat,tflat,trflat/100.,tflat/100.,fshift,efshift
ap1=ap1-fshift
ap=ap-fshift

if not keyword_set(bin) then begin
  ap1=ap1/3.*8.
  ap=rebin(ap,27,4096)
  delta1=delta1/3.*8.
endif

;force ap to be exactly as ap1 !!!
for i=0,n_elements(ap1)-1 do ap[i,*]=ap1[i]

;save aperture info
if keyword_set(bin) then mode='8x2' else mode='1x1'
writefits,strcompress('ap_'+mode+'.fits',/rem),ap
writefits,strcompress('ap1_'+mode+'.fits',/rem),ap1
openw,1,strcompress('delta1_'+mode+'.dat',/rem)
printf,1,delta1
close,1
printf,10,'delta1=',delta1
print,'delta1=',delta1



hbias, f
dispdir,f,idisp 
gain = sxpar(header,'GAIN')
rdnoise = sxpar(header,'RDNOISE')
f = f * gain
vf =  f + rdnoise^2
if keyword_set(bin) then begin
	f = rebin(f,514*nboost,2048)
	vf = rebin (vf, 514*nboost, 2048)
endif

collapse, f, idisp, ap, delta1, xf, vf= vf, vs=xfv
ws,'xflat.fits',xf,xfv, hd=header


printf,10,'fshift/efshift=',fshift,efshift
print,'fshift/efshift=',fshift,efshift


;remove scattered light and extract spes
printf,10,'removing scattered light and extracting spes ...'
print,'removing scattered light and extracting spes ...'
for i=0,n_elements(wspe)-1 do begin
  j=wspe[i]
  frame = readfits(st[j].filename,header)
  hbias,frame,rdn=rdn,bin=bin
  if keyword_set(bin) then delta2 = delta1/3. else delta2=delta1
  scatter,frame,idisp,delta2,back
  frame = frame - back
  gain=sxpar(header,'GAIN')
  rdnoise=sxpar(header,'RDNOISE')
  printf,10,st[j].filename,' rdnoise=',rdnoise/gain,rdn,' (counts)'
  print,st[j].filename,' rdnoise=',rdnoise/gain,rdn,' (counts)'
  frame = frame * gain
  vframe = frame + rdnoise^2
  if keyword_set(bin) then begin
	frame = rebin (frame, 514*nboost, 2048)
	vframe = rebin (vframe, 514*nboost, 2048)
  endif
  collapse, frame, idisp, ap, delta1, xframe, vf = vframe, vs = xvframe 
  ws, 'x'+st[j].filename, xframe, xvframe, hd=header
endfor

;extract cals
printf,10,'extracing cals ...'
print,'extracing cals ...'
for i=0,n_elements(wcal)-1 do begin
  j=wcal[i]
  frame = readfits(st[j].filename,header)
  hbias,frame,rdn=rdn,bin=bin
  gain=sxpar(header,'GAIN')
  rdnoise=sxpar(header,'RDNOISE')
  printf,10,st[j].filename,' rdnoise=',rdnoise/gain,rdn,' (counts)'
  print,st[j].filename,' rdnoise=',rdnoise/gain,rdn,' (counts)'
  frame = frame * gain
  vframe = frame + rdnoise^2
  if keyword_set(bin) then begin
	frame = rebin (frame, 514*nboost, 2048)
	vframe = rebin (vframe, 514*nboost, 2048)
  endif
  collapse, frame, idisp, ap, delta1, xframe, vf = vframe, vs = xvframe 
  ws, 'x'+st[j].filename, xframe, xvframe, hd=header
endfor


;wave cal. cals
printf,10,'wavelength calibration of cals ...'
print,'wavelength calibration of cals ...'
for i=0,n_elements(wcal)-1 do begin
  printf,10,'calibrating ... '+'x'+st[wcal[i]].filename
  print,'calibrating ... '+'x'+st[wcal[i]].filename
  rs, 'x'+st[wcal[i]].filename, xframe, xvframe, hd=header
  xcal, xframe, wframe, bin=bin
  help,xframe,xvframe,wframe
  ws, 'x'+st[wcal[i]].filename, xframe, xvframe, w = wframe, hd=header
  if i eq 0 then begin
    calfiles = 'x'+st[wcal[i]].filename
    calmjd = st[wcal[i]].mjd0
  endif else begin
    calfiles = [ calfiles, 'x'+st[wcal[i]].filename ]
    calmjd = [ calmjd, st[wcal[i]].mjd0 ]  
  endelse
endfor

;average flatfield (approx. cal. using last cal)
rs,'xflat.fits', xframe, xvframe, hd=header
ws, 'xflat.fits', xframe, xvframe, w = wframe, hd=header


print,calfiles
print,calmjd

;read instrumental response
res=readfits('response.fits',hd)

;spec
for i=0,n_elements(wspe)-1 do begin
  printf,10,'calibrating ... '
  print,'calibrating ... '
  printf,10,'x'+st[wspe[i]].filename+'  mjd=',st[wspe[i]].mjd0
  print,'x'+st[wspe[i]].filename+'  mjd=',st[wspe[i]].mjd0
  rs, 'x'+st[wspe[i]].filename, xframe, xvframe, norder=norder, hd=header
  creject, xframe, xframe2
  xframe = xframe2
  wpre = max(where(st[wspe[i]].mjd0-calmjd gt 0.))
  wpos = min(where(st[wspe[i]].mjd0-calmjd lt 0.))
  printf,10,'   1-'+calfiles[wpre]+'  mjd=',calmjd[wpre]
  printf,10,'   2-'+calfiles[wpos]+'  mjd=',calmjd[wpos]
  print,'   1-'+calfiles[wpre]+'  mjd=',calmjd[wpre]
  print,'   2-'+calfiles[wpos]+'  mjd=',calmjd[wpos]
  rs, calfiles[wpre], x1, xv1, w = w1
  rs, calfiles[wpos], x2, xv2, w = w2
  if max(wpre) gt -1 and max(wpos) gt -1 then begin
    d = (st[wspe[i]].mjd0 - calmjd[wpre] ) / (calmjd[wpos] - calmjd[wpre] )
  endif else begin
	if max(wpre) gt -1 and max(wpos) lt  0 then d=0.
 	if max(wpre) lt  0 and max(wpos) gt -1 then d=1.
	if max(wpre) lt  0 and max(wpos) lt  0 then d=-1.
  endelse
  if d lt 0. then begin
	printf,10,'ERROR -- cannot wave calibrate!' 
	print,'ERROR -- cannot wave calibrate!' 
  endif else begin
  	wframe = (1. - d) * w1 + d * w2 
  	ws, 'x'+st[wspe[i]].filename, xframe, xvframe, w = wframe, hd=header

	;normalize by the flatfield
	for j=0,norder-1 do begin
	  ratio = 1.0d0 ; median(xf[j,*])/median(xframe[j,*])
	  xframe[j,*] = xframe[j,*] / xf[j,*] * ratio 
	  xvframe[j,*] = xvframe[j,*] / xf[j,*]^2 * ratio^2
	endfor

	;normalize by response
        xframe = xframe / res
	xvframe = xvframe / res^2	

	;trim
	margin=240
	np=n_elements(xframe[0,*])
	xframe = xframe[*,margin:np-1-margin]
	xvframe = xvframe[*,margin:np-1-margin]
        wframe = wframe[*,margin:np-1-margin]

	;conti, xframe, con
	;xframe = xframe/con
	;xvframe = xvframe/con^2

	ws, 'n'+st[wspe[i]].filename, xframe, xvframe, w = wframe, hd=header

  endelse
endfor
close,1

end

