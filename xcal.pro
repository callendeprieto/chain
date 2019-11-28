pro  xcal,c,x,calstats=calstats,bin=bin

;
;	IN:  c -float array 	Collapsed ThAr spectrum
;	OUT: x -float array 	Wavelengths corresponding to flux array c
;
;	KEYWORDS: calstats   7x27 array with 7 columns giving: aperture, nlines, 
;					dispersion, yerror1, yerror2, yerror3, yerror4
;
;			where yerrorN is the uncertainty for a N-th order polynomial calib.
;			bin - integer	Set to indicate 8x2 binning


ptol=5.2 ; pixels from the expected location to accept a line fit
tol=0.07 ; nm around nominal central wavelength to search for a line in the HORuS spectrum
wtol=3.3 ; max fwhm accepted for good lines (units of pixels)

calstats=dblarr(7,27)
nx=4096
if keyword_set(bin) then nx=nx/2

x=dblarr(27,nx)
openr,11,'/home/callende/idl/hors/hors_thar.dat'
for i=0,26 do begin


  ;plot,c[i,*]

  i2=0
  nlines=0
  order=0
  readf,11,i2,nlines,order
  coef=dblarr(order+1)
  readf,11,coef
  xx=poly(dindgen(4088),coef)
  if keyword_set(bin) then xx=rebin(xx,2044)

  print,'aperture = ',i+1
  print,nlines,' candidate calibration lines'

  ;lets measure
  ngood=0
  pixels=-1.d0
  lambdas=-1.d0
  for j=0,nlines-1 do begin
   
     ;first McD atlas
     i2=0
     j2=0
     pix=0.0d0
     l0=0.0d0
     readf,11,i2,j2,pix,l0
     wx=where(abs(xx-l0) lt tol)

     ;now fit the HORuS data
     ;plot,dindgen(n_elements(wx)),c[i,wx]
     g=gaussfit(dindgen(n_elements(wx)),c[i,wx],nter=4,a)
     ;oplot,dindgen(n_elements(wx)),g,col=180


     nwx=n_elements(wx)
     ;help,nwx
     ;print,a
     ;print,nwx/2.-1.-ptol,nwx/2.-1.+ptol
     ;print,wtol

     if a[0] gt 100. and a[1] gt nwx/2.-1.-ptol and a[1] lt nwx/2.-1.+ptol and abs(a[2]) gt 1.0 and abs(a[2]) lt wtol then begin

     ;plot,findgen(nwx),c[i,wx],title=string(a[0])+string(a[1])+string(a[2]),psym=-2
     ;oplot,findgen(nwx),g,col=180
     ngood=ngood+1
     pixels=[pixels,a[1]+wx[0]]
     lambdas=[lambdas,l0]
 
    endif
    
  endfor
 
  print,ngood,' successfully measured lines'

  yerror1=0.0
  yerror2=0.0
  yerror3=0.0
  yerror4=0.0

  if ngood le 1 then begin

    calstats[*,i]=[i+1,0.,0.,$
	yerror1,yerror2,yerror3,yerror4]

    x[i,*]=dblarr(nx)

  endif else begin

    pixels=pixels[1:ngood]
    lambdas=lambdas[1:ngood]

    c1=poly_fit(pixels,lambdas,1,yerror=yerror1)
    if ngood gt 3 then c2=poly_fit(pixels,lambdas,2,yerror=yerror2)
    if ngood gt 4 then c3=poly_fit(pixels,lambdas,3,yerror=yerror3)  
    if ngood gt 5 then c4=poly_fit(pixels,lambdas,4,yerror=yerror4)
    print,i+1,n_elements(pixels),median(xx-shift(xx,1)),$
	  yerror1,yerror2,yerror3,yerror4

    if yerror4 gt 0.0 then xpixels=poly(pixels,c4) else begin
      if yerror3 gt 0.0 then xpixels=poly(pixels,c3) else begin
        if yerror2 gt 0.0 then xpixels=poly(pixels,c2) else begin
          xpixels=poly(pixels,c1)
        endelse
      endelse
    endelse

    yerror1=0.0
    yerror2=0.0
    yerror3=0.0
    yerror4=0.0

    wgood=where(abs(lambdas-xpixels) lt 0.003)
    ngood=n_elements(wgood)
    pixels=pixels[wgood]
    lambdas=lambdas[wgood]
    c1=poly_fit(pixels,lambdas,1,yerror=yerror1)
    if ngood gt 3 then c2=poly_fit(pixels,lambdas,2,yerror=yerror2)
    if ngood gt 4 then c3=poly_fit(pixels,lambdas,3,yerror=yerror3)  
    if ngood gt 5 then c4=poly_fit(pixels,lambdas,4,yerror=yerror4)
    print,i+1,n_elements(pixels),median(xx-shift(xx,1)),$
	yerror1,yerror2,yerror3,yerror4

    if yerror4 gt 0.0 then x[i,*]=poly(dindgen(nx),c4) else begin
      if yerror3 gt 0.0 then x[i,*]=poly(dindgen(nx),c3) else begin
        if yerror2 gt 0.0 then x[i,*]=poly(dindgen(nx),c2) else begin
          x[i,*]=poly(dindgen(nx),c1)
        endelse
      endelse
    endelse

    calstats[*,i]=[i+1,n_elements(pixels),median(x[i,*]-shift(x[i,*],1)),$
	yerror1,yerror2,yerror3,yerror4]

    if 1 eq 0 then begin
      if n_elements(pixels) gt 1 then begin
        !p.multi=[0,0,0]
        plot,pixels,lambdas,psy=-2
        plot,pixels,lambdas-poly(pixels,c4),yr=[-0.005,0.005],psy=2
        oplot,[0,1e4],[-0.002,-0.002],linestyle=2
        oplot,[0,1e4],[0.002,0.002],linestyle=2
        ;wait,3
      endif
    endif

  endelse

endfor
close,11

end

