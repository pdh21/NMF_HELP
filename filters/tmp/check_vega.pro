
; 'cam_lw.res'

   filt_read_vega, vega_lam, vega_flam
;approximate conversion into ergs s-1 cm-2 A-1
   vega_flam = vega_flam*2.0e-17


   vega_nu = 2.998e18/vega_lam  ; Hz  (Note c in A s-1)
   vega_fnu = vega_flam*vega_lam^2/2.998e18 ;fnu in ergs s-1 cm-2 Hz-1  (Note c in As-1 and lam in A)

   vega_fnu = vega_fnu*1.e-3    ; W m-2 Hz-1

   vega_nufnu = vega_nu*vega_fnu

   plot, vega_lam, vega_nufnu/1.e-8, xrange=[3000, 1.1e4], xs=1

;http://www.stsci.edu/cgi-bin/NICMOS/si.pl?nav=tools:tools_prop&sel=id:291

;http://archive.stsci.edu/copernicus/vega.html
   plot, vega_lam, vega_flam/vega_lam, xrange=[2000, 3200]
;
;http://www.jach.hawaii.edu/JACpublic/UKIRT/software/cgs4/cgs4dr/node22.html

   print, interpol(vega_flam, vega_lam, [0.5556, 7.8, 20.]*1.e4)*1.e4 ; ergs s-1 cm-2 mum-1
   print, interpol(vega_fnu, vega_lam, [0.5556, 7.8, 20.]*1.e4)

; Fig. 13.3 p. 397 Peacock

