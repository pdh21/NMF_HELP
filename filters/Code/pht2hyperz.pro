; converts phot filters into hyperz format


dir = '~/Research/iso/calibration/Filter/'

cd, dir

files = findfile('*_r.dat')

get_lun, unit
openw, unit, 'pht.res'

FOR i=0, n_elements(files)-1 DO BEGIN 

readcol, files[i], lam, trans
lam = lam*1.e4

   nmax = n_elements(lam)
   name = strupcase(strmid(files[i], 0, strpos(files[i], 'fm_r.dat')))

   printf, unit, nmax, 'ISO PHOT '+name, format='(i9,5x,a)'

   FOR j=0, nmax-1 DO begin
      printf,unit, j+1, lam[j], trans[j], $
       format='(i9,2x,f10.2,2x,f10)'
   ENDFOR
   
ENDFOR

free_lun, unit
END

