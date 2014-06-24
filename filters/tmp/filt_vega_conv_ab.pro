PRO filt_vega_conv_ab, filt_lam, filt_trans, convab, zero, $
                       vega_lam=vega_lam, vega_flam=vega_flam, surf=surf

; calculates the conversion between vega and AB mags for
; this filter
;
; convab = mAB - mVega = mAB(Vega) = -2.5 log Int R fVega dnu  -48.60 
;                                             ----------------
;                                             Int R  dnu
; from
; http://webast.ast.obs-mip.fr/hyperz/hyperz_manual1/node19.html#vegased
;
;
; The "Zero Point" for converting between Magnitudes and Jy 
; is then defined to be MAB(vega) - 8.9
;
;  
;
; 

   IF n_params() LT 2 THEN message, 'Calling Sequence:'+$
    ''

   IF NOT keyword_set(vega_lam) OR NOT  keyword_set(vega_flam) THEN BEGIN
      message, 'vega_lam or vega_flam not set reading in VEGA file', /inf
      filt_read_vega, vega_lam, vega_flam
      vega_flam = vega_flam*2.0e-17

   ENDIF

;----------------------------------------------------------------------
; standard normalisation
;----------------------------------------------------------------------

   filt_trans0 = filt_trans/max(filt_trans)
  
;----------------------------------------------------------------------


   vega_fnu = vega_flam * vega_lam^2 / 2.998e18
   vega_nu =  2.998e18/vega_lam

; surface function 

   IF NOT keyword_set(surf) THEN surf = int_trap(filt_lam, filt_trans0)

   filt_nu =  2.998e18/filt_lam
   surf_nu = int_trap(filt_nu, filt_trans0)

;
; estimating transmission at points where vega has flux
; 

   vega_trans = interpol(filt_trans0, filt_lam, vega_lam)

;
; estimating vega flux at points where transmission is defined
; 

   filt_vega_fnu = interpol(vega_fnu, vega_lam, filt_lam)

; setting to zero outside range for which filter is defined

   filt_lam_minmax = minmax(filt_lam)
   out = where(vega_lam LT filt_lam_minmax[0] OR $
               vega_lam GT filt_lam_minmax[1],  nout)
   in = where(vega_lam ge filt_lam_minmax[0] and $
               vega_lam le filt_lam_minmax[1],  nvega_in)


   IF nout GT 0 THEN vega_trans[out] = 0.

; calculating correction factor

; N.B. first estimate using f lambda was wrong because it was
; assuming an unconventional spectral slope 
; 
;
;   convab1 = -2.5 * alog10( int_trap(vega_lam, vega_trans*vega_flam) / surf) ; - 48.60

;
; there is some choice about the interation steps and the best choice
; will depend on the variation of the spectrum and the filter
; two simple choices are (1) using the steps defined in the reference
; spectrum or (2) the steps defined in the filter
; we decide between these on the basis of the number of points
; 

   IF nvega_in lT n_elements(filt_lam) THEN begin
; filter has more points so is assumed to be better defined
      convab = -2.5 * alog10( int_trap(filt_nu, filt_trans0*filt_vega_fnu) / surf_nu)  - 48.60
   ENDIF else BEGIN                                 
; vega has more points so is assumed to be better defined
      convab = -2.5 * alog10( int_trap(vega_nu, vega_trans*vega_fnu)      / surf_nu) - 48.60
   endelse

;print, convab2

   IF !debug EQ 1 THEN BEGIN 
      stop
      plot,vega_lam,vega_fnu,xrange=filt_lam_minmax 
      oplot, filt_lam, filt_vega_fnu, color=2
      plot,vega_lam,vega_trans*vega_fnu,xrange=filt_lam_minmax 
      oplot, filt_lam, filt_trans0*filt_vega_fnu, color=2

   ENDIF

   zero = convab-8.9
      
   return

END



