PRO FILT_read_next, unit, npoint, name, lam, trans,  silent=silent


; reads the next filter transmission from a FILTER.RES file (used in 
; e.g. HYPERZ)

   npoint = 0
   name = ''
   lam = 0.
   trans = 0.
   
   CASE  eof(unit) OF 

      1: BEGIN
         message, 'End of File', /inf
         close, unit
         unit = -1
         return
      ENDCASE
      
      0:BEGIN
         readf, unit, npoint, name
         lam = fltarr(npoint)
         trans = fltarr(npoint)
; next step might be better without a loop
         FOR i=0, npoint-1 DO BEGIN
            readf, unit, j, lam0, trans0
            lam[i] = lam0
            trans[i] = trans0 > 0.  ;N.B. revised to ensure transmission allways >= 0.
         ENDFOR
         
      ENDCASE
      ELSE: message, 'Should not be here'
      
   ENDCASE

; check that file is in increasing order of lam

   check_order = sort(lam)
   check = where(check_order - indgen(n_elements(lam)) NE 0,  nbad)
   IF nbad GT 0 THEN BEGIN
      get_lun, unit2
      openw, unit2, strcompress(name, /rem)
      lam = lam[check_order]
      trans = trans[check_order]
      FOR j=0, npoint-1 DO BEGIN
         printf,unit2, j+1, lam[j], trans[j], format='(i9,2x,f10.2,2x,f10)'
      ENDFOR
      free_lun, unit2
      message, 'RES file has wavelengths out of order: output to '+strcompress(name, /rem)   , /inf
   ENDIF
   
END



