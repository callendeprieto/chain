pro	dispdir,f,idisp

;
;	Identify  dispersion direction. 
;
;	IN:      f - float array	Input image
;	OUT: idisp - integer            Dispersion direction (1 or 2)
;
;

sf=size(f)
np=(sf[1]+sf[2])/2.


mf=mean(f)*np
f1=total(f,1)/mf
f2=total(f,2)/mf
std=[stddev(f1),stddev(f2)]
idisp=(where(std eq max(std))+1)[0]

end

