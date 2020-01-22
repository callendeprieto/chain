pro xwindows1,ap1,delta1,left1,right1

;determine extraction windows when the central apertures (ap1) are not trazed


nap=n_elements(ap[*,0])
np=n_elements(ap[0,*])

left1 = fix(ap1 - delta1)
right1 = fix(ap1 + delta1)

end
