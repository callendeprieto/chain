pro ws,file,s,vs,w=w,hd=hd
;+
;	Writes a FITS echelle spectrum
;
;	IN: 	file	- string	full name of the input .ec.fits file
;	
;		s	- double	flux matrix	(npix x norder)
;		vs      - double        variance matrix (npix x norder)
;
;       KEYWORD:  w       - double        wavelength matrix (npix x norder) 
;                 hd     - strarr       header
;
;-

if N_params() LT 1 then begin
      print,'% ws: - ws,infile,s,vs[,w=w]'
      return
endif

ss=size(s)
if ss[0] ne 2 then begin
	message,'spectrum must be 2D (npix x naper)'
	return
endif

if n_elements(w) gt 0 then begin
	f=dblarr(3,ss[2],ss[1])
	f[0,*,*]=transpose(w)
	f[1,*,*]=transpose(s)
	f[2,*,*]=transpose(vs)
endif else begin
	f=fltarr(2,ss[2],ss[1])
        f[0,*,*]=transpose(s)
        f[1,*,*]=transpose(vs)
endelse

writefits,file,f,hd

end

