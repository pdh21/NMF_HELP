readfmt, '~/idl/filters/wfc/Readme', '(x,i3,32x,a11)', skip=8, id, name

get_lun, unit
openw, unit, '~/idl/filters/wfc/filter_wfc.res'
plot, findgen(10), findgen(10), xrange=[100, 1000], yrange=[0, 100.], /nodata
FOR i=0, n_elements(id)-1 DO BEGIN 
   file = '~/idl/filters/wfc/f'+strcompress(id[i], /rem)
   readcol, file, lam, trans
   oplot, lam, trans

   order = sort(lam)
   lam = lam(order)
   trans = trans(order)
; ensuring +ve ness
   trans = trans > 0
   printf, unit, n_elements(lam), 'INT WFC '+name[i], format='(i9,5x,a)'
   FOR j=0, n_elements(lam)-1 DO printf, unit, j+1, lam[j]*10., trans[j], format='(i9,2x,f10.2,2x,f10)'

ENDFOR
free_lun, unit
END
