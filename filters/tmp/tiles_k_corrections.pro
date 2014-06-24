get_lun, unit
   openw, unit, 'tiles_ks.txt'


   cww_dir = '/starpc07_1/starpc07_data/sjo/hyperz1.1/ZPHOT/templates/'
   cww_file = ['CWW_E_ext', 'CWW_Sbc_ext', 'CWW_Scd_ext', 'CWW_Im_ext']


   FOR j=0, 0 DO BEGIN 

      readcol, cww_dir+cww_file[j]+'.sed', sed_lam2, sed_flam2, /silent

      dir = '~/idl/filters/data/'
      file_stem = 'filter_wfc'
      filt_file = dir+file_stem+'.res'


      FOR i=1, 21 DO BEGIN 
         FILT_read, file_stem, dir=dir, i, npoint, name, filt_lam, filt_trans, silent=1
         plot,filt_lam,filt_trans,title=name


         k_correction, sed_lam2, sed_flam2, filt_lam, filt_trans, 1.5, k_z
         filt_vega_mag, sed_lam2, sed_flam2, filt_lam, filt_trans,  mag
         printf,unit, format='(a25,a20,4f7.2)', name, cww_file[j], 2.5*alog10(k_z), mag, mag+ 2.5*alog10(k_z), mag+ 2.5*alog10(k_z)+19.5-9.11
         

      ENDFOR


      dir = '~/idl/filters/data/'
      file_stem = 'filter_ukirt'
      filt_file = dir+file_stem+'.res'


      FOR i=1, 6 DO BEGIN 
         FILT_read, file_stem, dir=dir, i, npoint, name, filt_lam, filt_trans, silent=1
         plot,filt_lam,filt_trans,title=name

         k_correction, sed_lam2, sed_flam2, filt_lam, filt_trans, 1.5, k_z
         filt_vega_mag, sed_lam2, sed_flam2, filt_lam, filt_trans,  mag
         printf,unit, format='(a25,a20,4f7.2)', name, cww_file[j], 2.5*alog10(k_z), mag, mag+ 2.5*alog10(k_z), mag+ 2.5*alog10(k_z)+19.5-9.11
         

      ENDFOR





   ENDFOR

free_lun, unit


END
