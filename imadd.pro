pro imadd,list,output,average=average

;
; adds a bunch of fits images whose file names are in the
; input string array list, and writes the result to disk
; 
;


for i=0,n_elements(list)-1 do begin
	d=readfits(list[i],hd)*1.d0
	if i eq 0 then tot=d else tot=tot+d
endfor
if keyword_set(average) then tot=tot/double(n_elements(list))
writefits,output,tot,hd


end

