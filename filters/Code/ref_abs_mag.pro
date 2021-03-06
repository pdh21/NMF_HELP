PRO ref_abs_mag, sed_lam2, sed_flam2, abs_mag, z, ref_filt, filts, omega=omega, h0=h0

   dir = '~/idl/filters/data/'
   file_stem = 'filter_seb'
   filt_file = dir+file_stem+'.res'

   IF NOT keyword_set(omega) THEN omega = 1
   IF NOT keyword_set(h0) THEN h0 = 100.
   IF NOT keyword_set(lambda) THEN lambda = 0




   FILT_read, file_stem, dir=dir, ref_filt, npoint, name, filt_lam, filt_trans, silent=1
   k_correction, sed_lam2, sed_flam2, filt_lam, filt_trans, z, k_z
   k_z = 2.5*alog10(k_z)
   filt_vega_mag, sed_lam2, sed_flam2, filt_lam, filt_trans,  mag_sed, $
    vega_lam=vega_lam, vega_flam=vega_flam
   mag_sed = mag_sed+k_z

   d_l = lum_dist(z, omega, h0=h0)
;N.B. conversion of d_l from  Mpc to 10pc
   mag_true = abs_mag+5*(alog10(d_l)+5.) + k_z
   
   ref_zero = mag_true-mag_sed

   print, 'Ref:',  name, 'M:', abs_mag, 'z:', z, 'm:', mag_true, $
    format='(a,2x,a,2x,a,2x,f7.2,2x,a,2x,f5.2,2x,a,2x,f6.2)'

   print, 'Filter', 'Vega', 'AB', 'muJy', format='(a25,2x,a5,2x,a5,2x,a5)'
   FOR i=0, n_elements(filts)-1 DO BEGIN 
      

      FILT_read, file_stem, dir=dir, filts[i], npoint, name, filt_lam, filt_trans, silent=1
      k_correction, sed_lam2, sed_flam2, filt_lam, filt_trans, z, k_z
      k_z = 2.5*alog10(k_z)
      filt_vega_mag, sed_lam2, sed_flam2, filt_lam, filt_trans,  mag_sed, $
       vega_lam=vega_lam, vega_flam=vega_flam
      mag_sed = mag_sed+k_z
      m = mag_sed+ref_zero
      
      filt_properties,filt_lam, filt_trans, $
       eff_lam=eff_lam, surf=surf, band=band, /silent
     
      filt_vega_conv_ab, filt_lam, filt_trans, convab, jy_zero, $
       vega_lam=vega_lam, vega_flam=vega_flam, surf=surf

;stop
;      convab = 0
;      zero = 0
;      print, convab, jy_zero
      print, name, m, m+convab, mag2flux(m, jy_zero-6*2.5), format='(a25,2x,f5.1,2x,f5.1,1x,f8.2)'

   ENDFOR

END

