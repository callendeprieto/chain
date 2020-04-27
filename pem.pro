pro	pem,y,n,yfit,percentile=percentile

;+
;	Break a spectrum in n pieces and normalize each 
;	by its own mean value.
;
;	IN: y	-fltarr		Input array
;		n	-integer	number of pieces 
;
;	OUT: yfit -fltarr	Output array
;
;	KEYWORDS: percentile	Desired percentile (0.5 for median)
;	
;	C. Allende Prieto, July 2011, T4 Barajas
;			 , January 2015, IAC -- added  percentile
;-

nel=n_elements(y)
np=nel/n
y2=y
yfit=y
for i=0,n-1 do begin
	if keyword_set(percentile) then mm=cmedian(y2[i*np:np*(i+1)-1],percentile=percentile) else $
	mm=mean(y2[i*np:np*(i+1)-1])
	yfit[i*np:np*(i+1)-1]=mm
	;y2[i*np:np*(i+1)-1]=y2[i*np:np*(i+1)-1]/mm
endfor
if np*i-1 lt nel-1 then begin
	if keyword_set(percentile) then mm=cmedian(y2[np*i:nel-1],percentile=percentile) else $
	mm=mean(y2[np*i:nel-1])
	yfit[np*i:nel-1]=mm
	;y2[np*i:nel-1]=y2[np*i:nel-1]/mm
endif

end
