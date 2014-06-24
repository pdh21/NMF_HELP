omega = 1
h0 = 100
filts = [253, 262, 260, 261, 246, 265, 266, 267, 268, 198, 199, 200, 201, 202, 203, 204]


;z = 1


cww_dir = '/starpc07_1/starpc07_data/sjo/hyperz1.1/ZPHOT/templates/'
cww_file = ['CWW_E_ext', 'CWW_Sbc_ext', 'CWW_Scd_ext', 'CWW_Im_ext']


j = 3
readcol, cww_dir+cww_file[j]+'.sed', sed_lam2, sed_flam2, /silent

; M* u
ref_filt = 164
abs_mag =-18.34
print, cww_file[j]
ref_abs_mag, sed_lam2, sed_flam2, abs_mag, z, ref_filt, filts, omega=omega, h0=h0


j = 1
readcol, cww_dir+cww_file[j]+'.sed', sed_lam2, sed_flam2, /silent

; M* g

ref_filt = 165
abs_mag =-20.4
print, cww_file[j]
ref_abs_mag, sed_lam2, sed_flam2, abs_mag, z, ref_filt, filts, omega=omega, h0=h0


j = 0
readcol, cww_dir+cww_file[j]+'.sed', sed_lam2, sed_flam2, /silent


; M* r
ref_filt = 166
abs_mag =-20.83
print, cww_file[j]
ref_abs_mag, sed_lam2, sed_flam2, abs_mag, z, ref_filt, filts, omega=omega, h0=h0


; M* i

ref_filt = 167
abs_mag =-21.26
print, cww_file[j]
ref_abs_mag, sed_lam2, sed_flam2, abs_mag, z, ref_filt, filts, omega=omega, h0=h0


; M* z

ref_filt = 168
abs_mag =-21.66
print, cww_file[j]
ref_abs_mag, sed_lam2, sed_flam2, abs_mag, z, ref_filt, filts, omega=omega, h0=h0
