;***********************************************************************
PRO filt_mag, sed_lam, sed_flam, filt_lam, filt_trans0,  mag

;+
; NAME:
;	filt_mag
;
; PURPOSE:
;	Calculates the magnitude for a given filter & SED
;
;
; CATEGORY:
;	analysis
;
; CALLING SEQUENCE:
;
;	filt_mag, sed_lam, sed_flam, filt_lam, filt_trans0, z, k_z
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
;IDL> filt_mag, sed_lam2, sed_flam2, filt_lam, filt_trans, z, k_z
;IDL> plot, z, 2.5*alog10(k_z), lines=0, thick=1
;
; MODIFICATION HISTORY:
; 	Written by:	Seb Oliver (Sussex) 24th Sep 2002
;	July, 1996	
;-


IF n_params() LT 5 THEN message, 'Calling sequence: filt_mag, sed_lam, sed_flam, filt_lam, filt_trans0, mag'


; setting up common wavelength axis

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

mag = -2.5*alog10(numerator[0])


END



