
;----------------------------------------------------------------------
; reading in filter
;----------------------------------------------------------------------

dir = '~/idl/filters/data/'
file_stem = 'filter_seb'
filt_file = dir+file_stem+'.res'


FILT_read, file_stem, dir=dir, 167, npoint, name, filt_lam, filt_trans, /silent

get_lun, unit
openw, unit, 'filter_167.txt'


;----------------------------------------------------------------------
; CWW SEDs
;----------------------------------------------------------------------

cww_dir = '/starpc07_1/starpc07_data/sjo/hyperz1.1/ZPHOT/templates/'
cww_file = ['CWW_E_ext',  'CWW_Sbc_ext', 'CWW_Scd_ext', 'CWW_Im_ext']

readfmt,'./andrew_barb.txt','(5x,a6,3x,f5)', id, z, skip=29, /silent


FOR i=0, n_elements(cww_file )-1 DO begin

   readcol, cww_dir+cww_file[i]+'.sed', sed_lam2, sed_flam2, /silent
   k_correction, sed_lam2, sed_flam2, filt_lam, filt_trans, z, k_z


   printf, unit, '----------------------------------------------------------------------'
   printf, unit, 'SPECTRUM: '+ cww_file[i]
   printf, unit, 'FILTER:   '+ name
   printf, unit, '----------------------------------------------------------------------'
   FOR j=0, n_elements(z)-1 DO printf, unit, id[j], z[j], 2.5*alog10(k_z[j])
   printf, unit, ''

ENDFOR
free_lun, unit


END
