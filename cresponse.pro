pro cresponse,hotstar,flat,response,margin=margin

;creates a response array from an extracted hot star spectrum
;first it normalizes the hot star spectrum using a flatfield
;second it runs a median filter on the ratio to create a correction 
;spectrum

macho=machar(/double)
tiny=macho.eps
buffer=100
margin=fix(2*buffer*1.2)

rs,hotstar,yh,vh
rs,flat,yf,vf

sf=size(yf)
if sf[0] ne 2 then begin
	message,'Dimension of data array must be 2!'
	return
endif


npix=sf[2]
norder=sf[1]
response=fltarr(norder,npix)

for i=0,norder-1 do begin
	a1=yf[i,*]
	a2=yh[i,*]
	response[i,*]=cmedian(a2/a1,buffer)
	response[i,0:margin]=1.0d0
	response[i,npix-margin:npix-1]=1.0d0
	plot,a2/a1,yr=[-2,2]
	oplot,response[i,*],col=180
	stop
endfor

end

