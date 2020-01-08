pro hfind,f,idisp,ap1,delta1,$
fwhm=fwhm,width=width,smoothinglength=smoothinglength,xorder=xorder,edge=edge

;
; Find approximately equidistant orders
;
;	IN:      f - float array		Input image
;		idisp - integer (1/2)		Dispersion direction
;	OUT: 	ap1 - float array		Locations of centers of the orders
;						in the cross-dispersion direction
;		delta1 - float                   Half width of the orders in the xdisp dir
;						(extraction should proceed between the
;						centers in the ap array and -/+ delta1)
;
;	Keywords: fwhm	- float			approx. FWHM of the orders 
;		  width - float			ratio between the width of the spectral 
;						orders and the distance between them 
;						(0.7 is a good choice for HORuS)
;		  smoothinglength -integer	length over which one should smooth the 
;						flux to help the peak search
;                 xorder - integer		order for the polynomial to fit the spacing
;						of the orders in the cross-disp. direction
;		  edge  - float			margin to the edges of the detector in the
;						cross-disp. direction to ignore orders 
;						in units of the average inter-order 
;						separation
;

sf=size(f)
np=sf[idisp] ; npixels in the dispersion direction
xdisp=1      
if idisp eq 1 then xdisp=2
nx=sf[xdisp] ; npixels in the cross-disp. direction


if not keyword_set(fwhm) then fwhm=10.
if not keyword_set(width) then width=0.7
if not keyword_set(smoothinglength) then smoothinglength=50
if not keyword_set(xorder) then xorder=4
if not keyword_set(edge) then edge=20.
margin=0.3 ; exclude margin*100 % of the peaks on each end to fit a polynomial to the locations of the orders in the x-disp. direction
margin=0.1
x=csmooth(total(f,idisp),smoothinglength)
;plot,x
peaks,x,loci,vals,fwhm=fwhm

nloci=n_elements(loci)
xo=findgen(nloci)
c=poly_fit(xo[nloci*margin:nloci*(1.-margin)],loci[nloci*margin:nloci*(1.0-margin)],xorder)
;plot,xo,loci,psy=4
;oplot,xo[nloci*margin:nloci*(1.0-margin)],loci[nloci*margin:nloci*(1.0-margin)],psy=4,symsize=4
;oplot,xo,poly(xo,c),col=180,thick=2,psy=2

print,'mean/std between the measured order locations and the fit=',mean(poly(xo[nloci*margin:nloci*(1.0-margin)],c)-loci[nloci*margin:nloci*(1.0-margin)]),stddev(poly(xo[nloci*margin:nloci*(1.0-margin)],c)-loci[nloci*margin:nloci*(1.0-margin)]),' pixels'

;print,xo,poly(xo,c)

w=where(poly(xo,c) gt 0. and poly(xo,c) lt nx) 
if max(w) lt 10 then begin
  print,'% FIND: -- error -- found less than 10 orders'
  return
endif

;remove extreme orders 
ap1=poly(xo[w],c)
delta1=median(ap1-shift(ap1,1))/2.*width
w=where(poly(xo,c) gt delta1*edge and poly(xo,c) lt nx-delta1*edge)
if max(w) lt 10 then begin
  print,'% FIND: -- error -- found less than 10 orders after removing extreme orders'
  return
endif


ap1=poly(xo[w],c)

print,'% FIND1: -- selected -- ',n_elements(w), ' orders' 

end
