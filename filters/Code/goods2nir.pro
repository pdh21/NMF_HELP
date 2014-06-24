; 200minutes 

;http://cfa-www.harvard.edu/cfa/oir/Research/irac/iracsens.html

;Frame        Point source
;Time     Sensitivity (5-sigma, in uJy)
;(sec)   3.6 um   4.5 um     5.8 um    8um
;-------------------------------------------
;200     3.1       5.0        16        23
;100     5.3       7.7        24        32
; 30      16        19        57        67
; 12      54        55       154       126
;  2     345       360      1040       450
;0.4    3850      4000     11500      4350


;3.6mu sensitivity for 200 minutes is 

 
; 0.4muJy


;Frame    Maximum unsaturated (80% of full well)
;Time           point source (milli-Jy) 
;(sec)   3.6 um    4.5 um    5.8 um     8um    Kmag
;-------------------------------------------    
;200      2.5        2.1        5.8      2.0    12.5
;100      5.1        4.2         12      4.5    11.7
; 30       17         14         39       16    10.4
; 12       43         35         98       41     9.4
;  2      256        212        590      250     7.5
;0.4     1280       1060       2950     1250     5.7
;0.02*  25000      21000      59000    25000     2.5

;saturation in 200s is 2.5mJy

;i.e. FOR a R-J star is 2.5 * (3.6/2.2)^2 =  6.7mJy 
;zero-point (Johnson K) =  -7.04
;
; == >  12.47

sat_36_mJy = [2.5, 5.1, 17., 43., 256., 1280, 25000]
sat_k_mJy = sat_36_mjy*(3.6/2.2)^2

print, flux2mag(sat_k_mjy*1.e-3, -7.04)


dbopen, '/starpc07_1/starpc07_data/sjo/2mass/N1_2mass'

k200 = dbfind('k_m<12.5')
k100 = dbfind('k_m<11.7', k200)
k30 = dbfind('k_m<10.4', k100)
k12 = dbfind('k_m<9.4', k30)
k2 = dbfind('k_m<7.5', k12)
k04 = dbfind('k_m<5.7', k2)
k002 = dbfind('k_m<2.5', k04)


dbext, k200, 'ra,dec', ra, dec
dbext, k100, 'ra,dec', ra1, dec1
dbext, k30, 'ra,dec', ra2, dec2
dbext, k12, 'ra,dec', ra3, dec3
dbext, k2, 'ra,dec', ra4, dec4
dbext, k04, 'ra,dec', ra5, dec5
;dbext, k002, 'ra,dec', ra6, dec6

survey_header, 'n1', n1, pfov=0.5



stringad, '16:09:40  +54:54:10', a, d

icmkhdr,  300, 300, 0.01, 0.01, 242.50500,  54.51000, 0, hd

im = readfits('/starpc07_1/starpc07_data/sjo/arch_maps/iras/n1/100.fit', hd)
hextract, im, hd, 100, 200, 100, 200

device, decom=0


set_plot, 'ps'
device, /col

colours
icplot, sigrange(im), hd, res=11, color=0, thick=5
icheader_oplot, hd, n1, thick=3, color=2, thick=5

icellipse_oplot, hd, a, d, 10./60., color=256
plotsym, 0, 0.6, /fill &  icellipse_oplot, hd, ra1*15., dec1, color=2, psym=8
plotsym, 0, 0.8, /fill & icellipse_oplot, hd, ra2*15., dec2, color=3, psym=8
plotsym, 0, 1, /fill & icellipse_oplot, hd, ra3*15., dec3, color=4, psym=8
plotsym, 0, 1.2, /fill & icellipse_oplot, hd, ra4*15., dec4, color=5, psym=8
plotsym, 0, 1,4 /fill & icellipse_oplot, hd, ra5*15., dec5, color=6, psym=8
;icellipse_oplot, hd, ra6*15., dec6, color=7
endps


$gv idl.ps & 

screen_grab,'tmp.bmp',/bmp 
