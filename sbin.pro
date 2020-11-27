pro sbin,im1,im2

;-
; Software binning 1x2 (4112x4096) -> 8x2 (514x2048)
;
; IN:   im1  -- string -- filename, unbinned image (1x1)
; OUT:  im2  -- string -- fiename, rebinned image  (8x2)
;
;+

d=readfits(im1,hd)
s=size(d)

if s[0] ne 2 or s[1] ne 4112 or s[2] ne 4096 then begin
  print,'SBIN: the input image ',im1,' does not have the expected dimentions'
  print,s
  return
endif

d2=rebin(d,s[1]/8,s[2]/2)
writefits,im2,d2,hd

end

