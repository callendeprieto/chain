pro chain,chaindir=chaindir,bin=bin,logfile=logfile

;basic dir/files
home=getenv('HOME')
if not keyword_set(chaindir) then chaindir=home+'/idl/chain'
if n_elements(bin) eq 0 then bin=1
if not keyword_set(logfile) then logfile='logfile'


if bin eq 1 then begin
  flatname='flat.fits'
  ap1name='ap1.fits'
  apname='ap.fits'
endif else begin
  flatname='uflat.fits'
  ap1name='uap1.fits'
  apname='uap.fits'
endelse

;params
scat='yes'  ; activate scattered light subtraction
la_minexptime=9.e10; min. exp. time (seconds) to activate la_cosmic removal
xmformat='vo'    ; format for order-merged output: '', 'vo', 'three-col' or 'rana'

;collect input files
files=file_search('00*fits')

;rebin unbinned data and rebin
if bin eq 0 then begin
  for i=0,n_elements(files)-1 do begin
    newfile=files[i]
    strput,newfile,'01',0
    sbin,files[i],newfile
  endfor
  files=file_search('01*fits')
endif

;make data inventory
inventory,files,st,/bin
if n_elements(st) eq 0 then begin
  if keyword_set(bin) then print,'% CHAIN: no binned (2x8) data' else $
	print,'% CHAIN: no unbinned (1x1) data'
  return
endif

openw,10,logfile

print,'make inventory of fits files ...'
print,'                      filename                object type  RA[deg] DEC[deg] t[s]'
print,'-------------------------------------------------------------------------------'
printf,10,'make inventory of fits files ...'
printf,10,'                      filename                object type RA[deg] DEC[deg] t[s]'
printf,10,'-------------------------------------------------------------------------------'
for i=0,n_elements(st.obstype)-1 do begin
        sfilename=strmid(st[i].filename,strpos(st[i].filename,'/',/reverse_search)+1,1000)
	printf,10,sfilename,st[i].object,strmid(st[i].obstype,0,3),$
		st[i].ra,st[i].dec,st[i].exptime,$
		format='(a42,x,a10,x,a3,x,f8.4,x,f7.3,x,i4)'
	print,sfilename,st[i].object,strmid(st[i].obstype,0,3),$
		st[i].ra,st[i].dec,st[i].exptime,$
		format='(a42,x,a10,x,a3,x,f8.4,x,f7.3,x,i4)'
endfor
print,'-------------------------------------------------------------------------------'
printf,10,'-------------------------------------------------------------------------------'

wcal=where(strmid(st.obstype,0,3) eq 'Cal')
wspe=where(strmid(st.obstype,0,3) eq 'Spe')
wfla=where(strmid(st.obstype,0,3) eq 'Fla')
wob=[wcal,wspe]


;have got enough data?
if max(wcal) lt 0 then begin
  print,'% CHAIN: no Cal images available, I quit' 
  return
endif
if max(wfla) lt 0 then begin
  print,'% CHAIN: no Fla images available, I quit' 
  return
endif
if max(wspe) lt 0 then begin
  print,'% CHAIN: no Spe images available, only Cal frames will collapsed' 
endif

;average binned flats and extract average
printf,10,'creating average binned flat ...'
print,'creating average flat ...'
imadd,st[wfla].filename,flatname,/av
f = readfits(flatname,header)

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
  endif
endif else begin
; or we take it from a reference image and use a ref flat to find the appropropriate offset
  ap1=readfits(chaindir+'/rap1.fits')
  ap=readfits(chaindir+'/rap.fits')

  ;derive spatial-direction offset of orders using flat
  rflat=readfits(chaindir+'/rflat.fits')
  trflat=total(rflat,2)
  tflat=total(f,2)
  ;if not keyword_set(bin) then tflat=rebin(tflat,514)
  xc,trflat,tflat,trflat/100.,tflat/100.,fshift,efshift
  ap1=ap1-fshift
  ;ap1=ap1/3.-fshift ;old ref. flat
  ap=ap-fshift
  width=0.72
endelse

;compute delta1 (order width) 
;it must be done consistently in the inspect/collapse (1/2) routines
delta1=median(ap1-shift(ap1,1))/2.*width
np=n_elements(ap[0,*])
delta=median(ap[*,np/2]-shift(ap[*,np/2],1))/2.*width



;save aperture info
writefits,strcompress(ap1name,/rem),ap1
writefits,apname,ap
printf,10,'delta1=',delta1
print,'delta1=',delta1
printf,10,'delta=',delta
print,'delta=',delta
printf,10,'fshift=',fshift
print,'fshift=',fshift

;find actual left/right limits for each aperture
xwindows,ap,delta,left,right

gain = sxpar(header,'GAIN')
rdnoise = sxpar(header,'RDNOISE')
f = f * gain
vf =  f + rdnoise^2

collapse, f, idisp, left,right, xf, vf= vf, vs=xfv ,/clean
ws,'x'+flatname,xf,xfv, hd=header


;remove cosmics, scattered light and extract spes
printf,10,'removing cosmics and scattered light and extracting spes ...'
print,'removing cosmics and scattered light and extracting spes ...'
for i=0,n_elements(wspe)-1 do begin
  print,i
  j=wspe[i]
  filename=st[j].filename
  frame = float(readfits(filename,header))
  if st[j].exptime gt la_minexptime then begin
    gain=sxpar(header,'gain')
    rdnoise=sxpar(header,'rdnoise')
    writefits,'tmp.fits',frame,header
    la_cosmic,'tmp.fits',gain=gain,readn=rdnoise,sigfrac=20.,sigclip=4.5
    frame = readfits('tmp-out.fits')
    file_delete,'tmp.fits','tmp-out.fits','tmp-mask.fits'
  endif
  hbias,frame,rdn=rdn,/bin
  if scat eq 'yes' then begin
    scatter,frame,idisp,delta1,back
    frame = frame - back
  endif
  gain=sxpar(header,'GAIN')
  rdnoise=sxpar(header,'RDNOISE')
  printf,10,filename,' rdnoise=',rdnoise/gain,rdn,' (counts)'
  print,filename,' rdnoise=',rdnoise/gain,rdn,' (counts)'
  frame = frame * gain
  vframe = frame + rdnoise^2
  collapse, frame, idisp, left, right, xframe, vf = vframe, vs = xvframe ,/clean
  xfilename='x'+strmid(filename,strpos(filename,'/',/reverse_search)+1)
  if strmid(xfilename,strlen(xfilename)-3,3) eq '.gz' then $
    xfilename=strmid(xfilename,0,strlen(xfilename)-3)
  ;upgrade header with HORuS coords (HORUSRA/HORUSDEC) and reduction time stamp
  hheader,header 
  ;write extracted file
  ws, xfilename, xframe, xvframe, hd=header
  
  ;check memory usage
  ;mem,used,tot,avail
  ;print,'memory used-total-avail (MB) = ',used,tot,avail
  ;printf,10,'memory used-total-avail (MB) = ',used,tot,avail
endfor

;extract cals
printf,10,'extracing cals ...'
print,'extracing cals ...'
for i=0,n_elements(wcal)-1 do begin
  j=wcal[i]
  filename=st[j].filename
  frame = readfits(filename,header)
  hbias,frame,rdn=rdn,/bin
  gain=sxpar(header,'GAIN')
  rdnoise=sxpar(header,'RDNOISE')
  printf,10,filename,' rdnoise=',rdnoise/gain,rdn,' (counts)'
  print,filename,' rdnoise=',rdnoise/gain,rdn,' (counts)'
  frame = frame * gain
  vframe = frame + rdnoise^2
  collapse, frame, idisp, left, right, xframe, vf = vframe, vs = xvframe ,/clean
  xfilename='x'+strmid(filename,strpos(filename,'/',/reverse_search)+1)
  if strmid(xfilename,strlen(xfilename)-3,3) eq '.gz' then $
       xfilename=strmid(xfilename,0,strlen(xfilename)-3)
  ;upgrade header with HORuS coords (HORUSRA/HORUSDEC) and reduction time stamp
  hheader,header 
  ;write extratcted file
  ws, xfilename, xframe, xvframe, hd=header

  ;check memory usage
  ;mem,used,tot,avail
  ;print,'memory used-total-avail (MB) = ',used,tot,avail
  ;printf,10,'memory used-total-avail (MB) = ',used,tot,avail
endfor


;wave cal. cals
printf,10,'wavelength calibration of cals ...'
print,'wavelength calibration of cals ...'
for i=0,n_elements(wcal)-1 do begin
  j=wcal[i]
  filename=st[j].filename
  xfilename='x'+strmid(filename,strpos(filename,'/',/reverse_search)+1)
  if strmid(xfilename,strlen(xfilename)-3,3) eq '.gz' then $
        xfilename=strmid(xfilename,0,strlen(xfilename)-3)
  printf,10,'calibrating ... '+xfilename
  print,'calibrating ... '+xfilename
  rs, xfilename, xframe, xvframe, hd=header
  xcal2, xframe, wframe,calstats=cs, chaindir=chaindir
  ratio=cs[6,*]/cs[2,*]
  minratio=min(ratio)
  maxratio=max(ratio)
  medratio=median(ratio)
  print,'min-max-median of rms/dispersion =',minratio,maxratio,medratio,format='(a30,3(1x,f10.5))'
  printf,10,'min-max-median of rms/dispersion =',minratio,maxratio,medratio,format='(a30,3(1x,f10.5))'

  ws, xfilename, xframe, xvframe, w = wframe, hd=header
  if i eq 0 then begin
    calfiles = xfilename
    calmjd = st[j].mjd0
  endif else begin
    calfiles = [ calfiles, xfilename ]
    calmjd = [ calmjd, st[j].mjd0 ]  
  endelse
endfor

;average flatfield (approx. cal. using last cal)
rs,'x'+flatname, xframe, xvframe, hd=header
ws, 'x'+flatname, xframe, xvframe, w = wframe, hd=header


print,calfiles
print,calmjd

;spec
for i=0,n_elements(wspe)-1 do begin
  j=wspe[i]
  filename=st[j].filename
  xfilename='x'+strmid(filename,strpos(filename,'/',/reverse_search)+1)
  if strmid(xfilename,strlen(xfilename)-3,3) eq '.gz' then $
       xfilename=strmid(xfilename,0,strlen(xfilename)-3)
  nfilename='n'+strmid(filename,strpos(filename,'/',/reverse_search)+1)
  if strmid(nfilename,strlen(nfilename)-3,3) eq '.gz' then $
       nfilename=strmid(nfilename,0,strlen(nfilename)-3)
  printf,10,'calibrating ... '
  print,'calibrating ... '
  printf,10,xfilename+'  mjd=',st[j].mjd0
  print,xfilename+'  mjd=',st[j].mjd0
  rs, xfilename, xframe, xvframe, norder=norder, hd=header
  creject, xframe, xframe2
  xframe = xframe2
  wpre = max(where(st[j].mjd0-calmjd gt 0.))
  wpos = min(where(st[j].mjd0-calmjd lt 0.))
  if max(wpre) gt -1 then begin  
    printf,10,'   1-'+calfiles[wpre]+'  mjd=',calmjd[wpre]
    print,'   1-'+calfiles[wpre]+'  mjd=',calmjd[wpre]
    rs, calfiles[wpre], x1, xv1, w = w1
    sxaddpar,header,'FILECAL1',calfiles[wpre],'Cal. spectrum preceding observation'
    sxaddpar,header,'MJDCAL1',calmjd[wpre],'MJD Ca. spectrum preceding'
  endif
  if max(wpos) gt -1 then begin  
    printf,10,'   2-'+calfiles[wpos]+'  mjd=',calmjd[wpos]
    print,'   2-'+calfiles[wpos]+'  mjd=',calmjd[wpos]
    rs, calfiles[wpos], x2, xv2, w = w2
    sxaddpar,header,'FILECAL2',calfiles[wpos],'Cal. spectrum following observation'
    sxaddpar,header,'MJDCAL2',calmjd[wpos],'MJD Ca. spectrum following'
  endif

  if max(wpre) gt -1 and max(wpos) gt -1 then begin
    d = (st[j].mjd0 - calmjd[wpre] ) / (calmjd[wpos] - calmjd[wpre] )
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
  	ws, xfilename, xframe, xvframe, w = wframe, hd=header

        ;normalize and order merge to create the m* files
        if xmformat ne '' then xm, xfilename,xmformat=xmformat

	;normalize by the flatfield
	for k=0,norder-1 do begin
	  ratio = 1.0d0 ; median(xf[k,*])/median(xframe[k,*])
	  xframe[k,*] = xframe[k,*] / xf[k,*] * ratio 
	  xvframe[k,*] = xvframe[k,*] / xf[k,*]^2 * ratio^2
	endfor

	;trim
	margin=150
	np=n_elements(xframe[0,*])
	xframe = xframe[*,margin:np-1-margin]
	xvframe = xvframe[*,margin:np-1-margin]
        wframe = wframe[*,margin:np-1-margin]

	ws, nfilename, xframe, xvframe, w = wframe, hd=header

  endelse
endfor
close,10

;create plots
;mkpng
cmd = 'python '+chaindir+'/horus.py'
spawn,cmd,back
print,back

end

