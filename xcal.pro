pro  xcal,c,x,calstats=calstats,chaindir=chaindir

;
;	IN:  c -float array 	Collapsed ThAr spectrum
;	OUT: x -float array 	Wavelengths corresponding to flux array c
;
;	KEYWORDS: calstats   7x27 array with 7 columns giving: aperture, nlines, 
;					dispersion, yerror1, yerror2, yerror3, yerror4
;
;			where yerrorN is the uncertainty for a N-th order polynomial calib.
;
;			chaindir - string  Path to the chain (both software and reference files)
;


ptol=8.0 ; pixels from the expected location to accept a line fit
tol=0.07 ; nm around nominal central wavelength to search for a line in the HORuS spectrum
wtol=3.3 ; max fwhm accepted for good lines (units of pixels)

calstats=dblarr(7,27)
nx=4096

home=getenv('HOME')
if not keyword_set(chaindir) then chaindir=home+'/idl/chain'

x=dblarr(27,nx)
;openr,11,chaindir+'/hors_thar.dat'
openr,11,chaindir+'/horus_thar.dat'
for i=0,26 do begin


  ;plot,c[i,*]

  i2=0
  nlines=0
  order=0
  readf,11,i2,nlines,order
  coef=dblarr(order+1)
  sig=dblarr(order+1)
  readf,11,coef
  readf,11,sig
  ;the following two lines are needed when using hors_thar.dat
  ;xx=poly(dindgen(4088),coef)
  ;xx=rebin(xx,2044)
  ;the following line needed when using horus_thar.dat
  xx=poly(dindgen(2048),coef)/10.

  print,'-----------------------------------------------------------'
  print,'aperture # ',i
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

     ;stop

     if a[0] gt 100. and a[1] gt nwx/2.-1.-ptol and a[1] lt nwx/2.-1.+ptol and abs(a[2]) gt 0.8 and abs(a[2]) lt wtol then begin

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

    calstats[*,i]=[i,0.,0.,$
	yerror1,yerror2,yerror3,yerror4]

    x[i,*]=dblarr(nx)

  endif else begin

    pixels=pixels[1:ngood]
    lambdas=lambdas[1:ngood]

    c1=poly_fit(pixels,lambdas,1,yerror=yerror1,sigma=sigma1)
    if ngood gt 3 then c2=poly_fit(pixels,lambdas,2,yerror=yerror2,sigma=sigma2)
    if ngood gt 4 then c3=poly_fit(pixels,lambdas,3,yerror=yerror3,sigma=sigma3)  
    if ngood gt 5 then c4=poly_fit(pixels,lambdas,4,yerror=yerror4,sigma=sigma4)
    print,'disp-rms (ord=1,2,3,4)=',median(xx-shift(xx,1)),$
	  yerror1,yerror2,yerror3,yerror4,format='(a24,5(1x,f9.5))'

    if yerror4 gt 0.0 then begin
      cb=c4
      s=sigma4
    endif else begin
      if yerror3 gt 0.0 then begin
        cb=c3 
        s=sigma3
      endif else begin
        if yerror2 gt 0.0 then begin
          cb=c2 
          s=sigma2
        endif else begin
          cb=c1
          s=sigma1
        endelse
      endelse
    endelse

    xpixels=poly(pixels,cb)

    yerror1=0.0
    yerror2=0.0
    yerror3=0.0
    yerror4=0.0

    wgood=where(abs(lambdas-xpixels) lt 0.009)
    ngood=n_elements(wgood)
    pixels=pixels[wgood]
    lambdas=lambdas[wgood]
    c1=poly_fit(pixels,lambdas,1,yerror=yerror1,sigma=sigma1)
    if ngood gt 3 then c2=poly_fit(pixels,lambdas,2,yerror=yerror2,sigma=sigma2)
    if ngood gt 4 then c3=poly_fit(pixels,lambdas,3,yerror=yerror3,sigma=sigma3)  
    if ngood gt 5 then c4=poly_fit(pixels,lambdas,4,yerror=yerror4,sigma=sigma4)
    print,'disp-rms (ord=1,2,3,4)=',median(xx-shift(xx,1)),$
	  yerror1,yerror2,yerror3,yerror4,format='(a24,5(1x,f9.5))'

    if yerror4 gt 0.0 then begin
      cb=c4
      s=sigma4
    endif else begin
      if yerror3 gt 0.0 then begin
        cb=c3 
        s=sigma3
      endif else begin
        if yerror2 gt 0.0 then begin
          cb=c2 
          s=sigma2
        endif else begin
          cb=c1
          s=sigma1
        endelse
      endelse
    endelse

    x[i,*]=poly(dindgen(nx),cb)

    if n_elements(cb) eq 5 then begin
      ca=transpose(cb)
      sa=s
    endif else begin
      ca=[transpose(cb),replicate(0.0,5-n_elements(cb))]
      sa=[s,replicate(1.0,5-n_elements(cb))]
    endelse


    if i eq 0 then begin
      cc=ca
      ss=sa 
    endif else begin
      cc=[[cc],[ca]]
      ss=[[ss],[sa]]
    endelse

    calstats[*,i]=[i,n_elements(pixels),median(x[i,*]-shift(x[i,*],1)),$
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

cc=transpose(cc)
ss=transpose(ss)

;enforce coherence in the calibration across orders
;calrefine,cc,ss

;for i=0,26 do begin
;  x[i,*]=poly(dindgen(nx),cc[i,*])
;endfor

end

