pro mkpng
;
;  loops over the chain output spectra to create 2d png plots
;
;

xn=file_search('[xn]*fits')
for i=0,n_elements(xn)-1 do begin
  set_plot,'Z'
  device,set_pixel_depth=24,decompose=2
  rs,xn[i],y,w=x,/plot,fcol=255,bcol=0
  write_png,strmid(xn[i],0,strlen(xn[i])-4)+'png',tvrd()
  device,/close
endfor
set_plot,'X'
end
