pro	inventory,st,verbose=verbose,bin=bin

;
;	This routine returns information on all the HORuS frames
;	present in the current dir 
;
;
;
;


files=file_search('0*fits')
if n_elements(files) eq 1 then begin
if strlen(files[0]) eq 0 then begin
  print,'% INVENTORY: Cannot find 0*fits files in the current directory'
  print,'%            Searching for data in ../*/0*fits  ...'
  files=file_search('../*/0*fits')
  if n_elements(files) eq 1 then begin
  if strlen(files[0]) eq 0 then begin
    print,'% INVENTORY: Cannot find 0*fits files'
    print,'%            Giving up  ...'
    return
  endif
  endif
endif
endif

for i=0,n_elements(files)-1 do begin

  info=file_info(files[i])

  ;print,info.size/1024./1024.

  if keyword_set(bin) then begin
;	if info.size/1024./1024. gt 2.1 then continue; MB
  endif else begin
	if info.size/1024./1024. lt 30. then continue; MB
  endelse

  ;read data and header
  print,'reading '+files[i]
  d=mrdfits(files[i],0,hd)

  ;extract info from header
  naxis=fix(sxpar(hd,'NAXIS'))
  naxis1=fix(sxpar(hd,'NAXIS1'))
  naxis2=fix(sxpar(hd,'NAXIS2'))
  bscale=double(sxpar(hd,'BSCALE'))
  bzero=double(sxpar(hd,'BZERO'))
  instrume=string(sxpar(hd,'INSTRUME'))
  object=string(sxpar(hd,'OBJECT'))
  exptime=double(sxpar(hd,'EXPTIME'))
  speed=float(sxpar(hd,'SPEED'))
  namps=fix(sxpar(hd,'NAMPS'))
  channel=string(sxpar(hd,'CHANNEL'))
  ccdsum=fix(sxpar(hd,'CCDSUM'))
  obstype=string(sxpar(hd,'OBSTYPE'))
  ra=double(sxpar(hd,'RADEG'))
  dec=double(sxpar(hd,'DECDEG'))
  rotang=double(sxpar(hd,'ROTANG'))
  ipa=double(sxpar(hd,'IPA'))
  lst=string(sxpar(hd,'LST'))
  mjd0=double(sxpar(hd,'MJD-OBS')) ; mjd at start of observation

  sti={filename: files[i], naxis: naxis, naxis1: naxis1, naxis2: naxis2, $
	bscale: bscale, bzero: bzero, $
        instrume: instrume, object: object, exptime: exptime, $
	speed: speed, namps: namps, channel: channel, ccdsum: ccdsum, $
	obstype: obstype, ra: ra, dec: dec, rotang: rotang, ipa: ipa, $
	lst: lst, mjd0: mjd0 }

  if n_elements(st) eq 0 then st=[sti] else st=[st,sti]

  ;stats
  ;d=d*bscale+bzero
  ;dmean=mean(d)
  ;dmedian=median(d)
  ;dstddev=stddev(d)

endfor


if n_elements(st) eq 0 then begin
  print,'% INVENTORY: -- no data available for the specified binning size'
  return
endif

;select only unbinned or binned data
if keyword_set(bin) then begin
        w=where(st.naxis1 eq 514 and st.naxis2 eq 2048)
        st=st[w]
endif else begin
        w=where(st.naxis1 eq 4112 and st.naxis2 eq 4096)
        st=st[w]
endelse


if keyword_set(verbose) then begin

	print,'make inventory of fits files ...'
	print,'                           filename                object  obstype    mjd0     exptime [s]'

	for i=0,n_elements(st.obstype)-1 do begin
  		print,st[i].filename,st[i].object,strmid(st[i].obstype,0,3),$
			st[i].mjd0,st[i].exptime,$
			format='(x,a45,x,a12,x,a3,x,f12.5,x,f7.0)'
	endfor
endif

end


