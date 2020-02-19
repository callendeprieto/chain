pro mkcube,list,d,layer=layer

;
; creates a 3D array from a list of 2D fits images
;
;  IN: list -- string	File with the names of the FITS images (one per line)
;
;  OUT: d -- fltarr     Output 3D data cube 
;
;  KEYWORD: layer -- when inputing extracted spectra (3D arrays)
;                      this parameter chooses the layer 
;			(for wavelength calibrated spectra 0 is wavelength,
;			1 is flux and 2 is variance; for uncalibrated spectra
;`			0 is flux and 1 is variance)
;

if n_elements(layer) eq 0 then layer=1

load,list,l

for i=0,n_elements(l)-1 do begin

    ;reading the data
    f=readfits(l[i],hd)
    sf=size(f)
    if sf[0] ne 2 and sf[0] ne 3 then begin
      print,'the size of the input arrays must be 2D (or 3D for extracted spectra)'
      print,sf
      return
    endif
    n1=sf[2]
    n2=sf[3]
    
    if i eq 0 then d=fltarr(n2,n1,n_elements(l))
    if sf[0] eq 2 then d[*,*,i]=f else d[*,*,i]=transpose(f[layer,*,*])

endfor

end


