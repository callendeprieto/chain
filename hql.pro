pro hql,im,x,y,palette=palette

;+
;  HORuS quick-look extraction
;
;
;-

if n_params() lt 1 then begin
  print,'% HQL: use -- hql,image[,palette=palette]'
  return
endif

if not keyword_set(palette) then palette=12

;read image, convert to e-, and compute variance
d=readfits(im,hd)
hbias,d
gain=sxpar(hd,'GAIN')
rdnoise=sxpar(hd,'RDNOISE')
d=d*gain
vd=d+rdnoise^2


n=sxpar(hd,'NAXIS')
n1=sxpar(hd,'NAXIS1')
n2=sxpar(hd,'NAXIS2')

if n ne 2 then begin
    print,'%HQG: -- ERROR -- the input HORuS image must be 2D (NAXIS=2)'
    return
endif

if n1 eq 514 and n2 eq 2048 then begin
  bin=1
endif else begin
  if n1 eq 514*8 and n2 eq 2048*2 then begin
    bin=0
  endif else begin
    print,'%HQG: -- ERROR -- the input HORuS image is neither in 1x1 nor in 8x2 binning'
    return
  endelse
endelse

;default apertures
idisp=2
ap1=readfits('/home/callende/idl/hors/rap1.fits',0,hd)
if bin eq 1 then ap1=ap1/3. else ap1=ap1/3./8.


;extract spectrum
collapse1,d,idisp,ap1, xd, vf=vd, vs=xvd

;read reference wavelength solution
norder=27
npix=4096
if bin then npix=npix/2
w=dblarr(norder,npix)
openr,11,'/home/callende/idl/hors/hors_thar.dat'
for i=0,26 do begin
  readf,11,i2,nlines,order
  coef=dblarr(order+1)
  readf,11,coef
  xx=poly(dindgen(4088),coef)
  xx=interpol(xx,npix)
  for j=0,nlines-1 do begin
    readf,11,i2,j2,pix,l0
  endfor
  w[i,*]=xx
endfor
close,11

order=1
step_1=median(w[*,0]-shift(w[*,0],1))
wobject=where(strpos(hd,'OBJECT') gt -1)
if min(wobject) gt -1 then wobject=min(wobject) else wobject=6

window,0,xsize=700,ysize=500
plot,findgen(npix),xd[order-1,*]/mean(xd[order-1,*])+order-1,$
    yrange=[0,norder+1],$
    ytitle='Aperture',title=hd(wobject),$
    charsize=1.2,xcharsize=.01,xrange=[-npix*0.2,npix*1.2],$
    xstyle=1,ystyle=1,col=0,back=255,position=[0.2,0.2,0.8,0.8]

loadct,palette
display,transpose(xd/sqrt(xvd)),/noerase,position=[0.15,0.1,0.8,0.8]
colorbar,col=10,charsi=1.5

xyouts,0.4,0.05,'Wavelength (nm)',/normal,charsize=1.5,col=0  

xyouts,-npix*0.28,order-1,w[order-1,0],charsi=1.3,col=0
xyouts,npix,order-1,w[order-1,npix-1],charsize=1.3,col=0
	
for order=1,norder do begin
  oplot,findgen(npix),xd[order-1,*]/mean(xd[order-1,*])+order-1,$
    col=0

  xyouts,-npix*0.28,order-1,w[order-1,0],charsi=1.3,col=0
  xyouts,npix,order-1,w[order-1,npix-1],charsize=1.3,col=0
endfor

x=w
y=xd

end
