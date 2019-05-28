pro getfwhm,xcal,file=file

root=strsplit(xcal,'.',/extract)
root=root[0]
if not keyword_set(file) then file=root+'.fwhm'

rs,xcal,ce,cev,w=x

loadct,12
u=20
openw,u,file
for i=0,n_elements(ce[*,0])-1 do begin
  plot,x[i,*],ce[i,*],chars=1.8,xtitle='Wavelength (nm)',ytitle='counts',title='ThAr'
  wait,1
  peaks,ce[i,*],loci,y,threshold=max(ce[i,*])/3.
  help,loci
  loci=loci[sort(y)]
  loci=loci[0:min([n_elements(loci)-1,20])]
  for j=0,n_elements(loci)-1 do begin
     if loci[j] lt 110 or loci[j] gt 4088-110 then continue
     plot,x[i,loci[j]-12:loci[j]+12],ce[i,loci[j]-12:loci[j]+12],psym=-2
     g=gaussfit(x[i,loci[j]-12:loci[j]+12],ce[i,loci[j]-12:loci[j]+12],nter=4,a)
     oplot,x[i,loci[j]-12:loci[j]+12],g,col=180
     a[2]=a[2]*2.0*sqrt(-2.0*alog(0.5))
     print,a
     printf,u,a
     ;wait,0.5
  endfor
endfor
close,u
end

