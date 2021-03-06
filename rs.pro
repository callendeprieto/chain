pro rs,infile,s,vs,w=w,norder=norder,npix=npix,hd=hd,plot=plot,fcol=fcol,bcol=bcol
;+
;	Reads a FITS echelle spectrum
;
;	IN: 	infile	- string	full name of the input .ec.fits file
;	
;	OUT:
;		s	- double	flux matrix	(npix x norder)
;               vs       - double        flux matrix     (npix x norder)
;
;	KEYWORDS: norder - integer	number of orders (xdisp dimension)
;		  npix	 - integer	number of pixels (dispersion dimen.)
;		  hd	 - strarr	header
;		  plot	 - 		produces a plot of the spectr. on demand
;                  w     - double        wavelength matrix (npix x norder)
;                  fcol  -  byte        color for lines in plot
;                  bcol  -  byte        color for background in plot
;
;-

if N_params() LT 1 then begin
      print,'% rs: - rs,infile,s,vs[,w=w,norder=norder,npix=npix,hd=hd,plot=plot]'
      return
endif

;colors
if n_elements(fcol) eq 0 then fcol=0
if n_elements(bcol) eq 0 then bcol=255

;reading the data
f=readfits(infile,hd)
sf=size(f)
if sf[0] ne 3 then begin
	message,'Dimension of data array must be 3!'
	return
endif

npix=sf[2]
norder=sf[3]

if n_elements(f[*,0,0]) gt 2 then begin
	w=transpose(f[0,*,*])
	s=transpose(f[1,*,*])
	vs=transpose(f[2,*,*])
endif else begin
        s=transpose(f[0,*,*])
        vs=transpose(f[1,*,*])
endelse

;picture it
if keyword_set(plot) then begin
	order=1
	;step_1=median(w[*,0]-shift(w[*,0],1))
	wobject=where(strpos(hd,'OBJECT') gt -1)
	if min(wobject) gt -1 then wobject=min(wobject) else wobject=6

	plot,findgen(npix),s[order-1,*]/mean(s[order-1,*])+order-1,$
	yrange=[0,norder+1],$
	ytitle='Aperture',title=hd(wobject),$
	charsize=1.2,xcharsize=.01,xrange=[-npix*0.2,npix*1.2],$
	xstyle=1,ystyle=1,col=fcol,back=bcol

        if n_elements(w) gt 0 then begin
	  xyouts,.42,.05,'Wavelength (nm)',/normal,charsize=1.5,col=fcol
	  xyouts,-npix*0.2,order,w[order-1,0],charsi=1.3,col=fcol
	  xyouts,npix,1,$
	  w[order-1,npix-1],charsize=1.3,col=fcol
        endif else begin
          xyouts,.45,.05,'Wavelength',/normal,charsize=1.5,col=fcol
        endelse

	
	for order=2,norder do begin
		oplot,findgen(npix),s[order-1,*]/mean(s[order-1,*])+order-1,$
			col=fcol
       		;xyouts,.4,.05,'Wavelength (Angstroms)',/normal,charsize=1.5,$
	        ;		col=fcol
                if n_elements(w) gt 0 then begin
        	  xyouts,-npix*0.2,order,w[order-1,0],charsi=1.3,col=fcol
        	  xyouts,npix,order,$
        	  w[order-1,npix-1],charsize=1.3,col=fcol
                endif

	endfor
endif

end

