pro snr,infile

;produces an image showing the S/N
;if the input is a extracted spectrum (x00*fits) it will read the data
;if the input is a raw one (00*fits) it will perform a quick extraction

home=getenv('HOME')
if not keyword_set(chaindir) then chaindir=home+'/idl/chain'

filename=strmid(infile,strpos(infile,'/',/reverse_search)+1)
if strmid(filename,0,3) eq 'x00' then begin
	rs,infile,s,vs
endif else begin      
        f = readfits(infile,header)  
	ap1=readfits(chaindir+'/rap1.fits')
        gain = sxpar(header,'GAIN')
        rdnoise = sxpar(header,'RDNOISE')
        f = f * gain
        vf =  f + rdnoise^2
        collapse1, f, 2, ap1-4.5,ap1+4.5, s, vf= vf, vs=vs
endelse

data=s/sqrt(vs)
loadct,12
display,transpose(data),charsize=1e-4
xyouts,0.02,0.07,'S/N',charsi=3,col=255,/norm
xyouts,0.05,0.02,infile,charsi=2,col=255,/norm
colorbar,maxrange=max(data),charsize=3,position=[0.12,0.16,0.9,0.24]
loadct,0

end
