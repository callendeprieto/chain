pro scombine,list,s

;
; combines a series of extracted spectra, removing cosmic rays and outliers
;
;


mkcube,list,d
load,list,l
first=readfits(l[0],hd)
rdnoise=sxpar(hd,'rdnoise')
gain=sxpar(hd,'gain')
cr_reject,d,rdnoise,0.0,1.0,gain,s,cnpix,cnoise,/noskyadjust

end

