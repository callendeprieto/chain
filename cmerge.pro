pro  cmerge,x,y,v,x2,y2,v2

;
; merge orders, adding the signal at overlapping wavelengths
;
;
;

sf=size(y)
if sf[0] ne 2 then begin
	message,'Dimension of data array must be 2!'
	return
endif

npix=sf[2]
norder=sf[1]

;deltay=2.*y[0,*]-shift(y[0,*],1)-shift(y[0,*],-1)
;wout=where(abs(deltay) gt 5.*stddev(deltay))
;if max(wout) gt -1 then y[0,wout]=(y[0,wout-1]+y[0,wout+1])*0.5
margin0=100
margin=margin0
x2=transpose(x[0,margin:npix-margin])
y2=transpose(y[0,margin:npix-margin])
v2=transpose(v[0,margin:npix-margin])
for i=1,norder-1 do begin
  margin=margin0-i*3 ; we cut more of the sides in the blue, reducing it towards the red
  xx=x[i,margin:npix-margin]
  yy=y[i,margin:npix-margin]
  vv=v[i,margin:npix-margin]
  ;deltay=2.*y[i,*]-shift(y[i,*],1)-shift(y[i,*],-1)
  ;wout=where(abs(deltay) gt 5.*stddev(deltay))
  ;if max(wout) gt -1 then y[i,wout]=(y[i,wout-1]+y[i,wout+1])*0.5
  ;plot,x[i-1,*],xr=[min(x[i-1,*]),max(x[i,*])],/xstyl,/nodata
  woverlap=where(xx le max(x2))
  if max(woverlap) gt -1 then begin
	wx=where(x2 ge min(xx))
	if max(xx) lt min(x2) then continue
        if max(wx) gt -1 then begin

                data=interpol(yy[woverlap],xx[woverlap],x2[wx])
                var=interpol(vv[woverlap],xx[woverlap],x2[wx])
                scale=mean(data/y2[wx])
                loadct,12
                plot,x2,y2,xr=[max(x2)-10.,max(xx[woverlap])+10.]
                oplot,x2[wx],data,col=100
                oplot,x2[wx],data/scale,col=180
                data=data/scale
                var=var/scale^2
		;y2[wx]=(y2[wx]/v2[wx]+interpol(y[i,woverlap]/v[i,woverlap],x[i,woverlap],x2[wx]))/$
		;(1./v2[wx]+interpol(1./v[i,woverlap],x[i,woverlap],x2[wx]))
		;y2[wx]=(y2[wx]+interpol(y[i,woverlap],x[i,woverlap],x2[wx]))/2.
		y2[wx]=(y2[wx]+data)/2.
                v2[wx]=v2[wx]+var/4.

		wnew=where(xx gt max(x2))
		;help,x2,x[i,wnew]

		x2=[x2,xx[wnew]]
		y2=[y2,yy[wnew]]
		v2=[v2,vv[wnew]]		
	endif
  endif else begin
	if max(x[i,*]) lt min(x2) then continue
	x2=[x2,transpose(xx)]
	y2=[y2,transpose(yy)]
	v2=[v2,transpose(vv)]
  endelse
  ;print,i,'---------'
  ;plot,x2,y2
endfor



end


