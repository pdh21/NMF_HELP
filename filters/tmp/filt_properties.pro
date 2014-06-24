PRO FILT_properties, lam, trans, $
                    eff_lam=eff_lam, surf=surf, band=band, silent=silent

; calculates some standard properties:
; based on the filter profile
    
;----------------------------------------------------------------------
; standard normalisation
;----------------------------------------------------------------------

   trans0 = trans/max(trans)
   
;----------------------------------------------------------------------
; effective lamlength
;----------------------------------------------------------------------

   eff_lam1 = int_trap(lam, trans0*lam)/int_trap(lam, trans0)


   IF NOT keyword_set(silent) THEN BEGIN
      eff_lam2 = int_tabulated(lam, trans0*lam)/int_trap(lam, trans0)
      eff_lam3 = int_tabulated2(lam, trans0*lam)/int_trap(lam, trans0)
      print, 'Effective Wavelength', eff_lam1, eff_lam2, eff_lam3
   ENDIF

   eff_lam = eff_lam1

;----------------------------------------------------------------------
; Surface function
;----------------------------------------------------------------------

   surf = int_trap(lam, trans0)
   IF NOT keyword_set(silent) THEN print, 'Surface function', surf

;----------------------------------------------------------------------
; Bandpass is the effective bandpass computed with a Gaussian approximation: 
;----------------------------------------------------------------------

   band_fn = (lam-eff_lam1)^2*trans0
   band = 2.*sqrt(int_trap(lam, band_fn)/surf)

   IF NOT keyword_set(silent) THEN print, 'Band pass',band
;----------------------------------------------------------------------
; conv_AB is the conversion between AB and Vega magnitudes (see Section A): 
;
;  Mvega=mab-conv_ab
;----------------------------------------------------------------------

END
