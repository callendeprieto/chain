pro chain1,logfile=logfile

if not keyword_set(logfile) then logfile='logfile'

;make data inventory 
inventory,st,/bin
if n_elements(st) eq 0 then begin
  print,'% CHAIN1: no binned (2x8) data'
  return
endif

openw,10,logfile


print,'make inventory of fits files ...'
print,'                               filename          object  obstype      RA        DEC          mjd0        exptime [s]'
printf,10,'make inventory of fits files ...'
printf,10,'                           filename                object  obstype      RA        DEC         mjd0        exptime [s]'
for i=0,n_elements(st.obstype)-1 do begin
	printf,10,st[i].filename,st[i].object,strmid(st[i].obstype,0,3),$
		st[i].ra,st[i].dec,st[i].mjd0,st[i].exptime,$
		format='(x,a45,x,a12,x,a3,x,f12.5,x,f12.5,x,f12.5,x,f7.0)'
	print,st[i].filename,st[i].object,strmid(st[i].obstype,0,3),$
		st[i].ra,st[i].dec,st[i].mjd0,st[i].exptime,$
		format='(x,a45,x,a12,x,a3,x,f12.5,x,f12.5,x,f12.5,x,f7.0)'
endfor


wcal=where(strmid(st.obstype,0,3) eq 'Cal')
wspe=where(strmid(st.obstype,0,3) eq 'Spe')
wfla=where(strmid(st.obstype,0,3) eq 'Fla')
wob=[wcal,wspe]

;have got enough data?
if max(wcal) lt 0 then begin
  print,'% CHAIN2: no Cal images available, I quit' 
  return
endif
if max(wfla) lt 0 then begin
  print,'% CHAIN2: no Fla images available, I quit' 
  return
endif
if max(wspe) lt 0 then begin
  print,'% CHAIN2: no Spe images available, only Cal frames will collapsed' 
endif

;average binned flats and extract average
printf,10,'creating average binned flat ...'
print,'creating average flat ...'
imadd,st[wfla].filename,'flat.fits',/av
f = readfits('flat.fits',header)

;find order information from flat...
printf,10,'find order information from flat spectrum...'
dispdir,f,idisp  
hbias,f,/bin

;we either find the apertures on the fly from the flat 
if (1 eq 0) then begin
  hfind,f,idisp,ap1,delta1,fwhm=4,width=0.65,smoothinglength=2,xorder=2
  if n_elements(ap1) ne 27 then begin
  	print,'% CHAIN: ERROR -- did not find exactly 27 orders!'
  	print,'               -- this is a requirement for xcal.'
  	printf,10,'% CHAIN: ERROR -- did not find exactly 27 orders!'
  	printf,10,'               -- this is a requirement for xcal.'
  	close,10
       stop
  endif
endif else begin
; or we take it from a reference image and use a ref flat to find the appropropriate offset
  ap1=readfits('/home/callende/idl/chain/rap1.fits')
  ap=readfits('/home/callende/idl/chain/rap.fits')

  ;derive spatial-direction offset of orders using flat
  rflat=readfits('/home/callende/idl/chain/rflat.fits')
  trflat=total(rflat,2)
  tflat=total(f,2)
  ;if not keyword_set(bin) then tflat=rebin(tflat,514)
  xc,trflat,tflat,trflat/100.,tflat/100.,fshift,efshift
  ap1=ap1-fshift
  ;ap1=ap1/3.-fshift ;old ref. flat
  width=0.72
endelse

;compute delta1 (order width) 
;it must be done consistently in the inspect/collapse (1/2) routines
delta1=median(ap1-shift(ap1,1))/2.*width



;save aperture info
writefits,strcompress('ap1.fits',/rem),ap1
printf,10,'delta1=',delta1
print,'delta1=',delta1
printf,10,'fshift=',fshift
print,'fshift=',fshift

;find actual left/right limits for each aperture
xwindows,ap1,delta1,left1,right1


gain = sxpar(header,'GAIN')
rdnoise = sxpar(header,'RDNOISE')
f = f * gain
vf =  f + rdnoise^2

help,f,idisp,left1,right1

collapse1, f, idisp, left1, right1, xf, vf= vf, vs=xfv
ws,'xflat.fits',xf,xfv, hd=header

;remove cosmics, scattered light and extract spes
printf,10,'removing cosmics and scattered light and extracting spes ...'
print,'removing cosmics and scattered light and extracting spes ...'
for i=0,n_elements(wspe)-1 do begin
  j=wspe[i]
  frame = readfits(st[j].filename,header)
  if st[j].exptime gt 120. then begin
    gain=sxpar(header,'gain')
    rdnoise=sxpar(header,'rdnoise')
    writefits,'tmp.fits',frame,header
    la_cosmic,'tmp.fits',gain=gain,readn=rdnoise
    frame = readfits('tmp-out.fits',header)
    file_delete,'tmp.fits','tmp-out.fits','tmp-mask.fits'
  endif
  hbias,frame,rdn=rdn,/bin
  scatter,frame,idisp,delta1,back
  frame = frame - back
  gain=sxpar(header,'GAIN')
  rdnoise=sxpar(header,'RDNOISE')
  printf,10,st[j].filename,' rdnoise=',rdnoise/gain,rdn,' (counts)'
  print,st[j].filename,' rdnoise=',rdnoise/gain,rdn,' (counts)'
  frame = frame * gain
  vframe = frame + rdnoise^2
  collapse1, frame, idisp, left1, right1, xframe, vf = vframe, vs = xvframe 
  ws, 'x'+st[j].filename, xframe, xvframe, hd=header
endfor

;extract cals
printf,10,'extracing cals ...'
print,'extracing cals ...'
for i=0,n_elements(wcal)-1 do begin
  j=wcal[i]
  frame = readfits(st[j].filename,header)
  hbias,frame,rdn=rdn,/bin
  gain=sxpar(header,'GAIN')
  rdnoise=sxpar(header,'RDNOISE')
  printf,10,st[j].filename,' rdnoise=',rdnoise/gain,rdn,' (counts)'
  print,st[j].filename,' rdnoise=',rdnoise/gain,rdn,' (counts)'
  frame = frame * gain
  vframe = frame + rdnoise^2
  collapse1, frame, idisp, left1, right1, xframe, vf = vframe, vs = xvframe 
  ws, 'x'+st[j].filename, xframe, xvframe, hd=header
endfor


;wave cal. cals
printf,10,'wavelength calibration of cals ...'
print,'wavelength calibration of cals ...'
for i=0,n_elements(wcal)-1 do begin
  printf,10,'calibrating ... '+'x'+st[wcal[i]].filename
  print,'calibrating ... '+'x'+st[wcal[i]].filename
  rs, 'x'+st[wcal[i]].filename, xframe, xvframe, hd=header
  xcal, xframe, wframe, /bin,calstats=cs
  print,'min(max) of rms/dispersion for 5th-order=',min(cs[6,*]/cs[2,*]),'(',$
                                      max(cs[6,*]/cs[2,*]),')'

  printf,10,'min(max) of rms/dispersion for 5th-order=',min(cs[6,*]/cs[2,*]),'(',$
                                      max(cs[6,*]/cs[2,*]),')'

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
  if max(wpre) gt -1 then begin  
    printf,10,'   1-'+calfiles[wpre]+'  mjd=',calmjd[wpre]
    print,'   1-'+calfiles[wpre]+'  mjd=',calmjd[wpre]
    rs, calfiles[wpre], x1, xv1, w = w1
  endif
  if max(wpos) gt -1 then begin  
    printf,10,'   2-'+calfiles[wpos]+'  mjd=',calmjd[wpos]
    print,'   2-'+calfiles[wpos]+'  mjd=',calmjd[wpos]
    rs, calfiles[wpos], x2, xv2, w = w2
  endif

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
        if abs(d) lt 1.e-7 then w2=w1
        if abs(d-1.) lt 1.e-7 then w1=w2
  	wframe = (1. - d) * w1 + d * w2 
  	ws, 'x'+st[wspe[i]].filename, xframe, xvframe, w = wframe, hd=header

	;normalize by the flatfield
	for j=0,norder-1 do begin
	  ratio = 1.0d0 ; median(xf[j,*])/median(xframe[j,*])
	  xframe[j,*] = xframe[j,*] / xf[j,*] * ratio 
	  xvframe[j,*] = xvframe[j,*] / xf[j,*]^2 * ratio^2
	endfor

	;trim
	margin=0
	np=n_elements(xframe[0,*])
	xframe = xframe[*,margin:np-1-margin]
	xvframe = xvframe[*,margin:np-1-margin]
        wframe = wframe[*,margin:np-1-margin]

	ws, 'n'+st[wspe[i]].filename, xframe, xvframe, w = wframe, hd=header

  endelse
endfor
close,10

end

