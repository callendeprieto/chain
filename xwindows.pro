pro xwindows,ap,delta,left,right

;determine extraction windows when the central apertures (ap) have been trazed


nap=n_elements(ap[*,0])

left=ap
right=ap

for i=0,nap-1 do begin
  left[i,*] = fix(ap[i,*] + delta) -fix(2.*delta)-1
  right[i,*] = fix(ap[i,*] + delta + 1)
endfor

end