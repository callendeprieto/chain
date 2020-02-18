pro mkrflat,sky,flat, chaindir=chaindir

;basic dir/files
home=getenv('HOME')
if not keyword_set(chaindir) then chaindir=home+'/idl/chain'


f = readfits(flat)
;o = readfits('0002326127-20191011-HORS-Spectroscopy.fits',hd)
s = readfits(sky,hd)

hbias,f,/bin
hbias,s,/bin

dispdir,f,idisp
ap1=readfits(chaindir+'/rap1.fits')
rflat=readfits(chaindir+'/rflat.fits')
trflat=total(rflat,2)
tflat=total(f,2)
xc,trflat,tflat,trflat/100.,tflat/100.,fshift,efshift
ap1=ap1/3.-fshift

;traze,f,idisp,ap1,ap,first=7,last=22,order=1
traze,s,idisp,ap1,ap,first=3,last=20,order=2,width=0.75

writefits,'rflat.fits',f
writefits,'rap1.fits',ap1
writefits,'rap.fits',ap

collapse,s,idisp,ap,s2
collapse1,s,idisp,ap1,s1
window,6
loadct,0
!p.multi=[0,1,2]
inspect,s[50:150,*],idisp,ap-50.
inspect1,s[50:150,*],idisp,ap1-50.,/noerase
window,1
loadct,0
!p.multi=[0,1,2]
inspect,s[200:350,*],idisp,ap-200.
inspect1,s[200:350,*],idisp,ap1-200.,/noerase
window,2
!p.multi=[0,1,2]
inspect,s[350:400,*],idisp,ap-350.
inspect1,s[350:400,*],idisp,ap1-350.,/noerase

window,6,xsize=1000,ysize=600
loadct,12
!p.multi=[0,0,0]
for i=0,26 do begin
  plot,s1[i,*],yr=[0,max(s2[i,*])*1.2]
  oplot,s2[i,*],col=180
  oplot,s1[i,*]
  wait,2
endfor
end
