pro mkpng
;
;  loops over the chain output spectra to create 2d png plots
;
;

xn=file_search('[xn]*fits')
dev = !D.Name
set_plot,'Z'
device,set_resolution=[512,256]*2, z_buffer=0
loadct,3
for i=0,n_elements(xn)-1 do begin
  rs,xn[i],y,w=x,/plot,fcol=0,bcol=200
  ;tvlct, r, g, b, /get
  write_png,strmid(xn[i],0,strlen(xn[i])-4)+'png',tvrd(/true)
  device,/close
endfor
set_plot,dev
end
