
;----------------------------------------------------------------------
; reading in filter
;----------------------------------------------------------------------

dir = '~/idl/filters/data/'
file_stem = 'filter_seb'
filt_file = dir+file_stem+'.res'

filts=[164,165,166,167,168]
filt_names=['u','g','r','i','z']

for ifilt=0,4 do begin 


FILT_read, file_stem, dir=dir, filts[ifilt], npoint, name, filt_lam, filt_trans, /silent

get_lun, unit
openw, unit, 'filter_'+filt_names[ifilt]+'.txt'


;----------------------------------------------------------------------
; CWW SEDs
;----------------------------------------------------------------------

cww_dir = '/starpc07_1/starpc07_data/sjo/Software/hyperz1.1/ZPHOT/templates/'
cww_dir = '/starpc07_1/starpc07_data/sjo/Software/hyperz1.1/ZPHOT/templates/'
cww_file = ['CWW_E_ext',  'CWW_Sbc_ext', 'CWW_Scd_ext', 'CWW_Im_ext']


;seds are defined from 1 AA thus they are valide for k-corrections up
; to z~8000 @  8140 AA the

z1 = 10.^(findgen(1000)/500.)

z = z1-1

FOR i=0, n_elements(cww_file )-1 DO begin

   readcol, cww_dir+cww_file[i]+'.sed', sed_lam2, sed_flam2, /silent
   k_correction, sed_lam2, sed_flam2, filt_lam, filt_trans, z, k_z


   printf, unit, '----------------------------------------------------------------------'
   printf, unit, 'SPECTRUM: '+ cww_file[i]
   printf, unit, 'FILTER:   '+ name
   printf, unit, '----------------------------------------------------------------------'
   FOR j=0, n_elements(z)-1 DO printf, unit, z[j], 2.5*alog10(k_z[j])
   printf, unit, ''

ENDFOR
free_lun, unit

endfor



END
