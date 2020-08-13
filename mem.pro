pro mem,used,tot,avail
;+
; Evaluates memory used by the current gdl/idl session, 
; the machine total, and the machine's available memory.
;
;  OUT: used - float  RAM memory used by gdl/idl (MB)
;       tot  - float  RAM memory in the computer (MB)  
;	used - float  RAM memory currently used (MB)
;-

help,/mem,output=memo
spawn,'free -b',free
memarr=strsplit(memo[0],' ',/ext)
freearr=strsplit(free[1],' ',/ext)

used=long64(strmid(memarr[3],0,strlen(memarr[3])-1))
tot=long64(freearr[1])
used2=long64(freearr[2])
avail=tot-used2

;from B to MB
used=used/1024.^2
tot=tot/1024.^2
avail=avail/1024.^2

end
