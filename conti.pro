pro conti,frame,con

;continuum normalize
sf=size(frame)
if sf[0] ne 2 then begin
	message,'Dimension of data array must be 2!'
	return
endif


npix=sf[2]
norder=sf[1]
con=dblarr(norder,npix)

for j=0,norder-1 do begin
	plot,frame[j,*]
        den = csmooth(transpose(cmedian(frame[j,*],40,perc=0.9)),200)
	;continuum,6,10,0.1,3.0,yc,den
	oplot,den,col=180,thick=2
	con[j,*]=den
	stop
endfor

end

