





get_lun, unit
openw, unit, 'filter_ukirt.res'

dir = '~/idl/filters/data/'

files = [$
'nsfcam_jmk_trans.dat', $
'nsfcam_hmk_trans.dat', $
'nsfcam_kmk_trans.dat', $
'nsfcam_kpmk_trans.dat', $
'nsfcam_lpmk_trans.dat', $
'nsfcam_mpmk_trans.dat']

names = [$
'UKIRT NSF CAM J ',$
'UKIRT NSF CAM H ',$
'UKIRT NSF CAM K ',$
'UKIRT NSF CAM K prime',$
'UKIRT NSF CAM L prime', $
'UKIRT NSF CAM M prime']

FOR j=0, n_elements(files)-1  DO BEGIN 

file = strcompress(files[j], /rem)


readcol, dir+file, lam, tran


lam = lam*1.e4

tran = tran > 0

order = sort(lam)

lam = lam[order]
tran = tran[order]


printf, unit, n_elements(lam), names[j], format='(i9,2x,a20)'
FOR i=0, n_elements(lam) -1 DO $
printf,unit, i+1, lam[i], tran[i], $ 
         format='(i9,2x,f10.2,2x,f10)'


ENDFOR
free_lun, unit



END


