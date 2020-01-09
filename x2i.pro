pro	x2i,xfile,ifile

;
; Converting an extracted HORuS spectrum, merging orders and resampling to a constant step.
; The output is readable in IRAF (splot).
;
;

rs,xfile,y,v,w=x,hd=hd
cmerge,x,y,v,x2,y2,v2

delta=5d-3
nfreq=(max(x2)-min(x2))/delta
xx=min(x2)+dindgen(nfreq)*delta

;help,x2,y2,xx
yy=interpol(y2,x2,xx)
sxaddpar,hd,'crpix1',1
sxaddpar,hd,'crval1',min(xx)*10.d0 ; AA
sxaddpar,hd,'cdelt1',delta*10.d0   ; AA
writefits,ifile,yy,hd

end
