;***********************************************************************
PRO filt_vega_mag, sed_lam, sed_flam, filt_lam, filt_trans0,  mag, $
                   vega_lam=vega_lam, vega_flam=vega_flam

;+
; NAME:
;	filt_vega_mag
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
;	filt_vega_mag, sed_lam, sed_flam, filt_lam, filt_trans0, mag, 
;                       [vega_lam=vega_lam, vega_flam=vega_flam]
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
;IDL> filt_vega_mag, sed_lam2, sed_flam2, filt_lam, filt_trans, z, k_z
;IDL> plot, z, 2.5*alog10(k_z), lines=0, thick=1
;
; MODIFICATION HISTORY:
; 	Written by:	Seb Oliver (Sussex) 19th Mar 2002
;	July, 1996	
;-






   IF NOT keyword_set(vega_lam) OR NOT  keyword_set(vega_flam) THEN BEGIN
      message, 'vega_lam or vega_flam not set reading in VEGA file', /inf
      filt_read_vega, vega_lam, vega_flam
      vega_flam = vega_flam*2.0e-17
 
   ENDIF


filt_mag, sed_lam, sed_flam, filt_lam, filt_trans0,  mag
filt_mag, vega_lam, vega_flam, filt_lam, filt_trans0,  vega_mag

mag = mag-vega_mag


END



