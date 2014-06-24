
;----------------------------------------------------------------------
; reading in filter
;----------------------------------------------------------------------

dir = '~/idl/filters/data/'
file_stem = 'filter_seb'
filt_file = dir+file_stem+'.res'

filts = [166, 164, 165, 166, 167,168, 136, 135, 137]

filts = [32, 12, 51, 32, 33,54, 34, 135, 35]

filts = [166, 32,  153]


   names = strarr(n_elements(filts))
   FOR j=0, n_elements(filts)-1 DO BEGIN 
      FILT_read, file_stem, dir=dir, filts[j], npoint, name, $
       filt_lam, filt_trans, /silent
      
      names[j] = strcompress(name)
   ENDFOR


print, '', names, format='(a10,15a12)'

;----------------------------------------------------------------------
; CWW SEDs
;----------------------------------------------------------------------

cww_dir = '/starpc07_1/starpc07_data/sjo/hyperz1.1/ZPHOT/templates/'
cww_dir = '~/ZPHOT/templates/'
cww_file = ['CWW_E_ext',  'CWW_Sbc_ext', 'CWW_Scd_ext', 'CWW_Im_ext']


FOR i=0, n_elements(cww_file )-1 DO begin

   readcol, cww_dir+cww_file[i]+'.sed', sed_lam, sed_flam, /silent

   mags = fltarr(n_elements(filts))

   FOR j=0, n_elements(filts)-1 DO BEGIN 
      FILT_read, file_stem, dir=dir, filts[j], npoint, name, $
      filt_lam, filt_trans, /silent


      filt_vega_mag, sed_lam, sed_flam, filt_lam, filt_trans,  mag, $
                   vega_lam=vega_lam, vega_flam=vega_flam

      IF keyword_set(ab) THEN BEGIN 
         filt_vega_conv_ab, filt_lam, filt_trans, convab, zero, $
          vega_lam=vega_lam, vega_flam=vega_flam, surf=surf
         mag = mag+convab
      endif

      mags[j] = mag
      
   ENDFOR
   mags = mags-mags[0]+20.

   print, cww_file[i], mags, format='(a10,15f12.2)'

ENDFOR


END
