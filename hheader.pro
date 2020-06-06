pro hheader,hd,voff=voff,chaindir=chaindir

;
; update header so that extracted spectra have info on the coordinates for the HORuS IFU
; (HORUSRA,HORUSDEC, both in degrees)
; and provide a timestamp and version info on the data reduction
;
; IN: hd -- strarr  Original header
;
; Keywords:
;     voff -- float  Velocity offset of the spectrum 

home=getenv('HOME')
if not keyword_set(chaindir) then chaindir=home+'/idl/chain'

hostname=getenv('HOST')
pwd=getenv('PWD')
user=getenv('USER')
version=''
openr,lun,chaindir+'/version.txt',/get_lun
readf,lun,version
close,lun
free_lun,lun

ra=sxpar(hd,'radeg')
dec=sxpar(hd,'decdeg')
mjd=sxpar(hd,'MJD-OBS')
jd=mjd+2400000.5d0

baryvel,jd,2000,vh,vb
vbary = vb[0]*cos(dec/!RADEG)*cos(ra*15.d0/!RADEG) + $   
   vb[1]*cos(dec/!RADEG)*sin(ra*15.d0/!RADEG) + vb[2]*sin(dec/!RADEG)


;HORuS IFU offset from GTC pointing (empirically derived, good to about 4 arcseconds)
hxoffset=-0.024576147  ; offset in ra*cos(delta) [deg], sigma=0.001 deg
hyoffset=-0.00392294   ; offset in dec [deg], sigma=0.001 deg


horusra=ra-hxoffset/cos(dec/180.d0*!dpi)
horusdec=dec-hyoffset

sxaddpar,hd,'horusra',horusra,'HORuS IFU RA deg'
sxaddpar,hd,'horusdec',horusdec,'HORuS IFU DEC deg'
sxaddpar,hd,'chain',version,'data reduction software version'
sxaddpar,hd,'chaindat',systime(0),'data reduction time stamp'
sxaddpar,hd,'chainman',user,'data reduction user'
sxaddpar,hd,'chaincom',hostname,'data reduction host'
sxaddpar,hd,'chainloc',pwd,'data reduction location'
sxaddpar,hd,'vbary',vbary,'projected velocity solar system barycenter'
if n_elements(voff) gt 0 then sxaddpar,hd,'vrad',voff+vbary

end
