get_lun, unit
   openw, unit, 'tiles_ks.txt'


   cww_dir = '/starpc07_1/starpc07_data/sjo/hyperz1.1/ZPHOT/templates/'
   cww_file = ['CWW_E_ext', 'CWW_Sbc_ext', 'CWW_Scd_ext', 'CWW_Im_ext']


   FOR j=0, 3 DO BEGIN 

      dir = '~/idl/filters/data/'
      file_stem = 'filter_wfc'
      filt_file = dir+file_stem+'.res'


      FOR i=1, 21 DO BEGIN 
         FILT_read, file_stem, dir=dir, i, npoint, name, filt_lam, filt_trans, silent=1
         plot,filt_lam,filt_trans,title=name
         readcol, cww_dir+cww_file[j]+'.sed', sed_lam2, sed_flam2, /silent

         k_correction, sed_lam2, sed_flam2, filt_lam, filt_trans, 1.5, k_z
         printf,unit, name, cww_file[j], 2.5*alog10(k_z)

      ENDFOR


      dir = '~/idl/filters/data/'
      file_stem = 'filter_ukirt'
      filt_file = dir+file_stem+'.res'


      FOR i=1, 6 DO BEGIN 
         FILT_read, file_stem, dir=dir, i, npoint, name, filt_lam, filt_trans, silent=1
         plot,filt_lam,filt_trans,title=name
         readcol, cww_dir+cww_file[j]+'.sed', sed_lam2, sed_flam2, /silent

         k_correction, sed_lam2, sed_flam2, filt_lam, filt_trans, 1.5, k_z
         printf,unit, name, cww_file[j], 2.5*alog10(k_z)

      ENDFOR





   ENDFOR

free_lun, unit


END
