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

x2=transpose(x[0,*])
y2=transpose(y[0,*])
v2=transpose(v[0,*])
for i=1,norder-1 do begin
  ;plot,x[i-1,*],xr=[min(x[i-1,*]),max(x[i,*])],/xstyl,/nodata
  woverlap=where(x[i,*] le max(x2))

  if max(woverlap) gt -1 then begin
	wx=where(x2 ge min(x[i,*]))
	if max(x[i,*]) lt min(x2) then continue
        if max(wx) gt -1 then begin

                data=interpol(y[i,woverlap],x[i,woverlap],x2[wx])
                var=interpol(v[i,woverlap],x[i,woverlap],x2[wx])
                scale=mean(data/y2[wx])
                loadct,12
                plot,x2,y2,xr=[max(x2)-10.,max(x[i,woverlap])+10.]
                oplot,x2[wx],data,col=100
                oplot,x2[wx],data/scale,col=180
                stop

                data=data/scale
                var=var/scale^2
		;y2[wx]=(y2[wx]/v2[wx]+interpol(y[i,woverlap]/v[i,woverlap],x[i,woverlap],x2[wx]))/$
		;(1./v2[wx]+interpol(1./v[i,woverlap],x[i,woverlap],x2[wx]))
		;y2[wx]=(y2[wx]+interpol(y[i,woverlap],x[i,woverlap],x2[wx]))/2.
		y2[wx]=(y2[wx]+data)/2.
                v2[wx]=v2[wx]+var/4.
		wnew=where(x[i,*] gt max(x2))
		;help,x2,x[i,wnew]
		x2=[x2,transpose(x[i,wnew])]
		y2=[y2,transpose(y[i,wnew])]
		v2=[v2,transpose(v[i,wnew])]		
	endif
  endif else begin
	if max(x[i,*]) lt min(x2) then continue
	x2=[x2,transpose(x[i,*])]
	y2=[y2,transpose(y[i,*])]
	v2=[v2,transpose(v[i,*])]
  endelse
  ;print,i,'---------'
  ;plot,x2,y2
endfor



end


