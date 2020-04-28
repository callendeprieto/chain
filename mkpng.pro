pro mkpng
;
;  loops over the chain output spectra to create 2d png plots
;
;

xn=file_search('[xn]*fits')
for i=0,n_elements(xn)-1 do begin
  set_plot,'z'
  device,decompose=0,set_resolution=[512,256]*2, z_buffer=0
  rs,xn[i],y,w=x,/plot,fcol=15,bcol=255
  write_png,strmid(xn[i],0,strlen(xn[i])-4)+'png',tvrd(/true)
  device,/close
endfor
set_plot,'X'
end
