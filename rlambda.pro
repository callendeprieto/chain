pro rlambda,np,x
;
;get a very rough wavelength calibration based on the orders' span
;  assuming linear dispersion
;
;  IN: np -- number of pixels in the disp. direction
;
load,'hors.txt',t,sk=7
t=float(t[*,0:26])
x=dblarr(27,np)
for i=0,26 do begin
  x[i,*]=interpol([t[1,i]-10.,t[3,i]+10.],np)-10.
endfor

end

