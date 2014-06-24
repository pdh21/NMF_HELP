;***********************************************************************
PRO k_correction, sed_lam, sed_flam, filt_lam, filt_trans0, z, k_z

;+
; NAME:
;	k_correction
;
; PURPOSE:
;	Calculates the K-correction for a given filter & SED
;
;
; CATEGORY:
;	analysis
;
; CALLING SEQUENCE:
;
;	k_correction, sed_lam, sed_flam, filt_lam, filt_trans0, z, k_z
;
; INPUTS:
;
;	 sed_lam: SED wavelength axis (in AA, does this matter)
;       sed_flam: SED flux density per unit wavelength (Flambda),
;                 arbitrary normalisation
;       filt_lam: Filter wavelength axis (in same units as sed_lam)
;       filt_trans0: Filter transmision (arbitrary normalisation)
;
; OPTIONAL INPUTS:
;
;              z: Array of redshifts
;	
;	
; KEYWORD PARAMETERS:
;	
;
; OUTPUTS:
;
; 	    k_z: K-correction as function of z
;
; OPTIONAL OUTPUTS:
;	
;              z: Array of redshifts
;
; COMMON BLOCKS:
;	
;
; SIDE EFFECTS:
;	
;
; RESTRICTIONS:
;	
;
; PROCEDURE:
;	
;
; impliments equation 13.14 (p.396) of Peacock.
;
; EXAMPLE:
;	
;IDL> dir = '~/idl/filters/data/'
;IDL> file_stem = 'filter_wfc'
;IDL> filt_file = dir+file_stem+'.res'
;IDL> 
;IDL> 
;IDL> FILT_read, file_stem, dir=dir, 1, npoint, name, filt_lam, filt_trans, silent=silent
;IDL> plot,filt_lam,filt_trans,title=name
;IDL> cww_dir = '/starpc07_1/starpc07_data/sjo/hyperz1.1/ZPHOT/templates/'
;IDL> cww_file = ['CWW_E_ext', 'CWW_Sbc_ext', 'CWW_Scd_ext', 'CWW_Im_ext']
;IDL> readcol, cww_dir+cww_file[2]+'.sed', sed_lam2, sed_flam2
;IDL> k_correction, sed_lam2, sed_flam2, filt_lam, filt_trans, z, k_z
;IDL> plot, z, 2.5*alog10(k_z), lines=0, thick=1
;
; MODIFICATION HISTORY:
; 	Written by:	Seb Oliver (Sussex) 19th Mar 2002
;	July, 1996	
;-




;----------------------------------------------------------------------
; impliments equation 13.14 (p.396) of Peacock.
;----------------------------------------------------------------------

IF n_params() LT 6 THEN message, 'Calling sequence: k_correction, sed_lam, sed_flam, filt_lam, filt_trans0, z, k_z'

IF n_elements(z) EQ 0 THEN BEGIN
   z = findgen(1000)/100.
   message, 'No input z, Default z set to 0<z<10, steps 0.01', /inf
ENDIF


;----------------------------------------------------------------------
; setting up common wavelength axis
; with fixed log lambda intervals 
; based on best of inputs
; Don't be tempted to use the actual SED wavelengths values within the 
; filter interval as the SED will be shifted to different redshifts Dooh!
;----------------------------------------------------------------------

; wavelength limits only need to be those of the filter

lambda_limits = alog(minmax(filt_lam))
lambda_range = lambda_limits[1]-lambda_limits[0]

; resolution needs to be the best of either

log_sed_lam = alog(sed_lam)
log_filt_lam = alog(filt_lam)
dlam_sed = ts_diff(log_sed_lam, 1)
dlam_filt = ts_diff(log_filt_lam, 1)

dlam_sed  = dlam_sed (0:(n_elements(dlam_sed)-2))
dlam_filt = dlam_filt(0:(n_elements(dlam_filt)-2))

dbin = min([abs(dlam_sed), abs(dlam_filt)])

; build up array


npoints = fix(lambda_range / dbin)+2

log_lam = findgen(npoints)*dbin+lambda_limits[0]

;----------------------------------------------------------------------
; estimating transmission at wavelength points
;----------------------------------------------------------------------

  log_filt_lam =  alog(filt_lam)
  filt_trans = interpol(filt_trans0,log_filt_lam, log_lam)

; enforce possitivity

  filt_trans = filt_trans > 0

;----------------------------------------------------------------------
; numerator
;----------------------------------------------------------------------

sed_lamflam0 = sed_lam*sed_flam

sed_lam_flam =  interpol(sed_lamflam0, log_sed_lam, log_lam)

numerator = int_trap(log_lam, filt_trans * sed_lam_flam)

;----------------------------------------------------------------------
; denominator
;----------------------------------------------------------------------


denom = fltarr(n_elements(z))
z1 = alog(1.+z)

for iz=0, n_elements(z)-1 DO BEGIN 
   sed_lam_flam =  interpol(sed_lamflam0, log_sed_lam, log_lam-z1[iz]) > 0

   denom[iz] = int_trap(log_lam, filt_trans * sed_lam_flam)

ENDFOR

k_z = numerator/denom

END



