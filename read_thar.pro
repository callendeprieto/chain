pro read_thar,wi,we,w,f,plot=plot,highest=highest,chaindir=chaindir
;+
;
;	Extracts a piece of the Th-Ar hollow cathode lamp used with the
;	2dcoude spectrograph (Allende Prieto 2001) 
;	
;	IN: wi - float  - Initial wavelength (angstroms)
;	    we - float  - Ending wavelength
;
;	OUT:w  - fltarr - Wavelength 
;	    f  - fltarr - Normalized flux
;	
;	KEYWORDS: plot	- produces a plot of the selected piece of the spectrum
;		  highest - labeling of the emission lines is restricted to 
;				5 strongest features
;                 chaindir - path to the ThAr FITS-format atlas data
;
;
;
;
;	C. Allende Prieto, UT, November 2001
;                       ", IAC, June 2020, changed data format from xdr to fits
;                                          modified for HORuS, added chaindir
;-
if N_params() LT 1 then begin		
      print,'% READ_THAR: Syntax - read_thar,lambda_i,lambda_e,wavelength,flux[,/plot]'
      return
endif
if (we gt 10596.401 or wi lt 3611.8806) then begin
	print,'% READ_THAR: the requested piece exceeds the boundaries'
	print,'% READ_THAR: of the atlas  3611.8806 -  10596.401 angstroms'
	return
endif
;restore,'~/idl/idl_database/thar/thar.xdr'

home=getenv('HOME')
if not keyword_set(chaindir) then chaindir=home+'/idl/chain'

data=readfits(chaindir+'/thar.fits')
wav=transpose(data[0,*])
flux=transpose(data[1,*])
mcd_thar_lines_file=chaindir+'/thar_output_combined.dat'
n=file_lines(mcd_thar_lines_file)
openr,lun,mcd_thar_lines_file,/get_lun
d=dblarr(3,n)
readf,lun,d
close,lun
free_lun,lun
f=flux[where(wav ge wi and wav le we)]
w=wav[where(wav ge wi and wav le we)]

if keyword_set(plot) then begin
	xb=min(w) & xr=max(w)
	if (max(where(d(0,*) gt xb and d(0,*) le xr)) gt -1) then begin
		dd=d[*,where(d(0,*) gt xb and d(0,*) le xr)]
		
		if keyword_set(highest) then begin
			noltww=n_elements(dd(0,*))
			maxlabels=5
			if (noltww gt maxlabels) then begin
				centralpix=where(min(abs(dd(1,0)-w)) eq abs(dd(1,0)-w))
				height=f[centralpix(0)]
				for iouy=1,noltww-1 do begin
					centralpix=where(min(abs(dd(1,iouy)-w)) eq abs(dd(1,iouy)-w))
					height=[height,f[centralpix[0]]]
				endfor
				dd2=dd
				pick=where(height eq max(height))
				height[pick(0)]=0.
				dd2=dd2[1,pick[0]]
				for iouy=1,maxlabels-1 do begin
					pick=[pick,where(height eq max(height))]
					height[pick[iouy]]=0.
				endfor
				dd=dd[*,pick]
			endif
		endif
		annot=1
	endif	
	
		; some parameters for plotting/annoting lines
	height=max([max(f)-min(f),500.])
	topband=1./2.2*height
	
	; plot
	
	plot,w,f,xtitle='Wavelength ('+string("305B)+')',thick=1,charsize=1.3,$
	xrange=[xb-1.,xr+1.0],xstyle=9,yr=[0,height+topband],ystyle=9,yminor=4
	
	; annotation
	if (annot ne 0) then begin
		length=1./3.*topband
		sep=1./3.*length
		extra=0.05 ; angstroms
		son=1.18 ; sizeofnumbers
		for i=0,n_elements(dd(0,*))-1 do begin
			
			; variable length of markers
			;heigh_low=min([max(y(where(x gt dd(1,i)-0.1 and $
			;x lt dd(1,i)+0.1)))+20*sep,height+sep])
			;oplot,[dd(1,i),dd(1,i)],[heigh_low,height+sep+length]
			
			; or constant length
			oplot,[dd(1,i),dd(1,i)],[height+sep,height+sep+length]
			
			if (dd(2,i) gt 0.01)  then begin
				xyouts,dd(1,i)+extra,height+1.5*sep+length,string(dd(1,i),$
				format='(f8.2)'),/data,charsize=son,orienta=90
			endif else begin
				if (dd(2,i) gt 0.001) then begin
					xyouts,dd(1,i)+extra,height+1.5*sep+length,string(dd(1,i),$
					format='(f9.3)'),/data,charsize=son,orienta=90
				endif else begin
					xyouts,dd(1,i)+extra,height+1.5*sep+length,string(dd(1,i),$
					format='(f10.4)'),/data,charsize=son,orienta=90
				endelse
			endelse
			
		endfor
	endif		
	
endif
end
