pro mem,used,tot,avail
;+
; Evaluates memory used by the current gdl/idl session, 
; the machine total, and the machine's available memory.
;
;  OUT: used - int  RAM memory used by gdl/idl (KB)
;       tot  - int  RAM memory in the computer (KB)  
;	used - int  RAM memory currently used (KB)
;-

help,/mem,output=memo
spawn,'free -b',free
memarr=strsplit(memo[0],' ',/ext)
freearr=strsplit(free[1],' ',/ext)

used=long(strmid(memarr[3],0,strlen(memarr[3])-1))
tot=long(freearr[1])
used2=long(freearr[2])
avail=tot-used2

end
