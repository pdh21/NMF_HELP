; converts a CAM CCGLWSPEC calibration file into the
; ascii format required for HYPERZ and other programmes.

no_zero = 1
cd, '/users/sjo/Research/iso/calibration/APR00'



cds_file = 'ccglwspec_98040514330660.cub'




restore, cds_file, /verb
help,cds,/str

help, data,/str


FOR i=0, n_elements(data)-1 DO print,data[i].EWHL, data[i].SWHL, data[i].PFOV, data[i].FCVF, data[i].WAVELENG, data[i].BANDWIDT

;000-055  8.063 --> 4.89  mum
;     308      88      60      95      6.75000      3.30000  LW2
;     308      88      60     110      11.5000      6.90000  LW10
;     308      88      60     125      15.0000      5.50000  LW3
;     308      88      60     140      15.0000      1.70000  LW9
;150-209  17.34 --> 11.58 mum CVF1 (159-209)
;210-235  11.48 --> 9.003 mum CVF2 (210-233)
;     308      88      60     245      11.4000      1.20000  LW8
;     308      88      60     260      9.63000      2.52000  LW7
;     308      88      60     275      7.75000      1.52000  LW6
;     308      88      60     290      6.75000     0.570000  LW5
;     308      88      60     305      6.00000      1.18000  LW4
;     308      88      60     320      4.50000      1.05000  LW1
;330-359 9.58-- > 8.116 mum  CVF?


name = 'ISO CAM '+['LW'+strcompress(indgen(10)+1, /rem), $
        'CVF0-'+string(indgen(56), format='(i3.3)'), $
        'CVF1-'+string(indgen(60)+150, format='(i3.3)'), $
        'CVF2-'+string(indgen(26)+210, format='(i3.3)'), $
        'CVF3-'+string(indgen(30)+330, format='(i3.3)')] 

fcvf = [320, 95, 125, 305,290, 275, 260, 245, 140, 110, $
        indgen(56), indgen(60)+150, indgen(26)+210, indgen(30)+330]

ident = intarr(n_elements(name))

FOR i=0, n_elements(name)-1 DO ident[i] = where(data.fcvf EQ fcvf[i])


start = 201

;----------------------------------------------------------------------
get_lun, unit
openw, unit, 'cam_lw.log'
FOR i=0, n_elements(name)-1 DO BEGIN
   len = strlen(name[i]) < 60
   format='(i6,1x,a'+strcompress(len, /rem)+','+strcompress(61-len, /rem)+'x,i6)'
  printf, unit, start+i, name[i], data[ident[i]].n, format=format
ENDFOR

;123456x123456789012345678901234567890123456789012345678901234567890x123456
;   193  M ISAAC (ESO web pages)                                        131   
free_lun, unit


;----------------------------------------------------------------------

;123456789xxxxx1234567890123456789012345678901234567890
;       13     Koo-Kron U+ filter (Koo s thesis)
;123456789xx1234567890xx123456789012345678901234567890
;        1     3000.00     0.00000

get_lun, unit
openw, unit, 'cam_lw.res'
FOR i=0, n_elements(name)-1 DO BEGIN
   nmax = data[ident[i]].n
   IF no_zero THEN temp = where(data[ident[i]].trans GT 0, nmax)
   printf, unit, nmax, name[i], format='(i9,5x,a)'
   FOR j=0, data[ident[i]].n-1 DO begin
      IF  data[ident[i]].trans[j] GT 0 THEN BEGIN 
         printf,unit, j+1, data[ident[i]].wave[j]*1.e4, data[ident[i]].trans[j], $
          format='(i9,2x,f10.2,2x,f10)'
      endif
   endfor
ENDFOR
free_lun, unit
END
