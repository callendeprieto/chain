pro creject,x,x2

; cosmic reject on extracted spectra

x2=x
nap=n_elements(x[*,0])
np=n_elements(x[0,*])
for i=0,nap-1 do begin
	m10=cmedian(x[i,*],10)
	m20=cmedian(x[i,*],20)
	sig=stddev(m10-m20)
	wcr=where(x[i,*] gt m20+10.*sig)
	if max(wcr) gt -1 then begin
		rest=where(x[i,*] le m20+10.*sig)
		x2[i,*]=interpol(x[i,rest],rest,dindgen(np))
	endif
	;plot,x[i,*]
	;oplot,m20,col=180
	;oplot,m20+10.*sig,col=100
	;oplot,wcr,x[i,wcr],psy=2
endfor

end

