pro mmatrix_FLS_z_lt_1
;TEST PROGRAM FOR RUNNING NMF_SPARSE ON PHOTOMETRY
 ;USES PARTS OF ISAAC'S SED FITTING CODE FOR DEALING WITH FILTERS
 
;specify the filters required (number corresponds to number in filter_seb.log)
 infilt=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,184,185,186,187,188,157,141,222,294,295,296]
 ;directory of filter_seb.log and filter_seb.res
;dir='/export/scratch/rh273/saves/'
filt='filter_master'


;Define redshift range and number of bins
tred=findgen(2)/2.0
nz=2.0
nif=n_elements(infilt)  ;number of infrared filters
fiwav=fltarr(nif,2501) ;lambda array for all infra-red filters
fitran=fltarr(nif,2501)  ;transmission array for all infrared filters


;load far-ir temps
restore, '/export/scratch/rh273/saves/mmodels.sav'
nit=n_elements(mmatrix[0,*]); number of templates
iwav=fltarr(nit,n_elements(lambdaGRID))
for i=0,nit-1 DO iwav[i,*]=lambdaGRID*10.0^4 ; lambda grid for each model, in angstroms
itflx=transpose(mmatrix)
for i=0,nit-1 DO itflx[i,*]=(itflx[i,*]*3*10.0^18)/iwav[i,*]^2 ; convert to SED flux density per unit wavelength (as required by filt_mag)

;load filters
irfeff=fltarr(nif) ; effective lambda for all infrared filters
irfsurf=fltarr(nif) ; surface function for all infrared filters
for k=0,n_elements(infilt)-1 do begin
;uses Seb's filt_read, loops for one res file, over filter ids
filt_read,filt,dir=dir,infilt[k],np,name,f1wav,f1tran
fiwav[k,0:np-1]=f1wav[0:np-1] 
fitran[k,0:np-1]=f1tran[0:np-1]/max(f1tran); normalised transmission
filt_properties,f1wav,f1tran,eff_lam=teff,surf=tsurf ;calculate effective lambda and surface function
irfeff[k]=teff
irfsurf[k]=tsurf
endfor

itmp=intarr(nit) ; number of infra red templates
irgrid=fltarr(nit,nif,nz) ; array number of templates by number of filters, by number of redshifts
 
 for i=0,nit-1 do begin
 for ii=0,nif-1 do begin
 for iii=0,nz-1 do begin
tz=tred[iii]
red=tred[iii]
tpts=where(iwav[i,*] ne 0,ntpts) ; where lambda template /= 0
;creates arrays that exclude values for when lambda =0
tiwav=fltarr(ntpts)
tiflx=fltarr(ntpts)
fpts=where(fiwav[ii,*] ne 0,nfpts)
tfwav=fltarr(nfpts)
tftran=fltarr(nfpts)
;fills arrays for values for when lambda /=0
tiwav[0:ntpts-1]=iwav[i,tpts]
tiflx[0:ntpts-1]=[itflx[i,tpts]]
tfwav[0:nfpts-1]=fiwav[ii,fpts]/(1.+tz) ;move filter to actual emission wavelength
tftran[0:nfpts-1]=fitran[ii,fpts]


filt_mag,tiwav,tiflx,tfwav,tftran,tmp;get magnitude for template i through filter ii

irgrid[i,ii,iii]=10^((tmp)*(-0.4))/irfsurf[ii]*(irfeff[ii]/(1.+red))^2/2.99792e18;convert flux(freq)
 
;tred[k]=tz
itmp[i]=i
;print,tred(k),cold_kcor(k),sb9_kcor(k)
endfor
endfor
ENDFOR
;save the irgrid file
save, irgrid,tred,nz, filename='~/idl/NMF/Rachels_Code_Bits/irgrid_test.sav'

;save, irgrid,tred,nz, filename='/export/scratch/rh273/saves/irgrid.sav'
END

pro prepare_data_FLS_z_lt_1
restore, '/export/scratch/rh273/saves/irgrid.sav'
;define number of photometric bands you are running NMF over
nbands=31
;load up data
dat=mrdfits('~/idl/subset-cosmos-april2012-xid-spire-match.fits',1)
;select redshift range

ind=where(dat.redshift gt 0.1 and dat.redshift lt 1.0,nsrc)
;dat=dat[ind]

;this line just checked data contained at least one SPIRE detection
;ind=where(dat.f250 gt 20.0 or dat.f350 gt 20.0 or dat.f500 gt 20.0, nsrc)


sed=fltarr(nsrc,nbands)
esed=fltarr(nsrc,nbands)

;Load up sed array with fluxes; fluxes in micro janskys

sed[*,0]=(3631/(10^((dat[ind].B)/2.5)))*10.0^6
no_data=where(dat[ind].B lt 0.0)
sed[no_data,0]=-99.0
sed[*,1]=(3631/(10^((dat[ind].V)/2.5)))*10.0^6
no_data=where(dat[ind].V lt 0.0)
sed[no_data,1]=-99.0
sed[*,2]=(3631/(10^((dat[ind].gp)/2.5)))*10.0^6
no_data=where(dat[ind].gp lt 0.0)
sed[no_data,2]=-99.0
sed[*,3]=(3631/(10^((dat[ind].rp)/2.5)))*10.0^6
no_data=where(dat[ind].rp lt 0.0)
sed[no_data,3]=-99.0
sed[*,4]=(3631/(10^((dat[ind].ip)/2.5)))*10.0^6
no_data=where(dat[ind].ip lt 0.0)
sed[no_data,4]=-99.0
sed[*,5]=(3631/(10^((dat[ind].zp)/2.5)))*10.0^6
no_data=where(dat[ind].zp lt 0.0)
sed[no_data,5]=-99.0
sed[*,6]=(3631/(10^((dat[ind].IB427)/2.5)))*10.0^6
no_data=where(dat[ind].IB427 lt 0.0)
sed[no_data,6]=-99.0
sed[*,7]=(3631/(10^((dat[ind].IB464)/2.5)))*10.0^6
no_data=where(dat[ind].IB464 lt 0.0)
sed[no_data,7]=-99.0
sed[*,8]=(3631/(10.0^((dat[ind].IB484)/2.5)))*10.0^6
no_data=where(dat[ind].IB484 lt 0.0)
sed[no_data,8]=-99.0
sed[*,9]=(3631/(10.0^((dat[ind].IB505)/2.5)))*10.0^6
no_data=where(dat[ind].IB505 lt 0.0)
sed[no_data,9]=-99.0
sed[*,10]=(3631/(10.0^((dat[ind].IB527)/2.5)))*10.0^6
no_data=where(dat[ind].IB527 lt 0.0)
sed[no_data,10]=-99.0
sed[*,11]=(3631/(10.0^((dat[ind].IB574)/2.5)))*10.0^6
no_data=where(dat[ind].IB574 lt 0.0)
sed[no_data,11]=-99.0
sed[*,12]=(3631/(10.0^((dat[ind].IB624)/2.5)))*10.0^6
no_data=where(dat[ind].IB624 lt 0.0)
sed[no_data,12]=-99.0
sed[*,13]=(3631/(10.0^((dat[ind].IB679)/2.5)))*10.0^6
no_data=where(dat[ind].IB679 lt 0.0)
sed[no_data,13]=-99.0
sed[*,14]=(3631/(10.0^((dat[ind].IB709)/2.5)))*10.0^6
no_data=where(dat[ind].IB709 lt 0.0)
sed[no_data,14]=-99.0
sed[*,15]=(3631/(10.0^((dat[ind].IB738)/2.5)))*10.0^6
no_data=where(dat[ind].IB738 lt 0.0)
sed[no_data,15]=-99.0
sed[*,16]=(3631/(10.0^((dat[ind].IB767)/2.5)))*10.0^6
no_data=where(dat[ind].IB767 lt 0.0)
sed[no_data,16]=-99.0
sed[*,17]=(3631/(10.0^((dat[ind].IB827)/2.5)))*10.0^6
no_data=where(dat[ind].IB827 lt 0.0)
sed[no_data,17]=-99.0
sed[*,18]=(3631/(10.0^((dat[ind].NB711)/2.5)))*10.0^6
no_data=where(dat[ind].NB711 lt 0.0)
sed[no_data,18]=-99.0
sed[*,19]=(3631/(10.0^((dat[ind].NB816)/2.5)))*10.0^6
no_data=where(dat[ind].NB816 lt 0.0)
sed[no_data,19]=-99.0
sed[*,20]=(3631/(10.0^((dat[ind].ks)/2.5)))*10.0^6
no_data=where(dat[ind].ks lt 0.0)
sed[no_data,20]=-99.0
sed[*,21]=(3631/(10.0^((dat[ind].u_s)/2.5)))*10.0^6
no_data=where(dat[ind].u_s lt 0.0)
sed[no_data,21]=-99.0
sed[*,22]=(3631/(10.0^((dat[ind].g_s)/2.5)))*10.0^6
no_data=where(dat[ind].g_s lt 0.0)
sed[no_data,22]=-99.0
sed[*,23]=(3631/(10.0^((dat[ind].r_s)/2.5)))*10.0^6
no_data=where(dat[ind].r_s lt 0.0)
sed[no_data,23]=-99.0
sed[*,24]=(3631/(10.0^((dat[ind].i_s)/2.5)))*10.0^6
no_data=where(dat[ind].i_s lt 0.0)
sed[no_data,24]=-99.0
sed[*,25]=(3631/(10.0^((dat[ind].z_s)/2.5)))*10.0^6
no_data=where(dat[ind].z_s lt 0.0)
sed[no_data,25]=-99.0
sed[*,26]=(3631/(10.0^((dat[ind].F814W)/2.5)))*10.0^6
no_data=where(dat[ind].F814W lt 0.0)
sed[no_data,26]=-99.0
sed[*,27]=(dat[ind].FLUX_24MIC)*10.0^3
no_data=where(dat[ind].FLUX_24MIC lt 0.0)
sed[no_data,27]=-99.0
sed[*,28]=(dat[ind].f250)*10.0^3
no_data=where(dat[ind].f250 lt 0.0)
sed[no_data,28]=-99.0
sed[*,29]=(dat[ind].f350)*10.0^3
no_data=where(dat[ind].f350 lt 0.0)
sed[no_data,29]=-99.0
sed[*,30]=(dat[ind].f500)*10.0^3
no_data=where(dat[ind].f500 lt 0.0)
sed[no_data,30]=-99.0

;Load up error sed array with fluxes; fluxes in micro janskys


esed[*,0]=(((alog(10))*(dat[ind].dB)*(dat[ind].B))/2.5)*10.0^6
no_data=where(dat[ind].dB lt 0.0)
esed[no_data,0]=-99.0
esed[*,1]=(((alog(10))*(dat[ind].dV)*(dat[ind].V))/2.5)*10.0^6
no_data=where(dat[ind].dV lt 0.0)
esed[no_data,1]=-99.0
esed[*,2]=(((alog(10))*(dat[ind].dgp)*(dat[ind].gp))/2.5)*10.0^6
no_data=where(dat[ind].dgp lt 0.0)
esed[no_data,2]=-99.0
esed[*,3]=(((alog(10))*(dat[ind].drp)*(dat[ind].rp))/2.5)*10.0^6
no_data=where(dat[ind].drp lt 0.0)
esed[no_data,3]=-99.0
esed[*,4]=(((alog(10))*(dat[ind].dip)*(dat[ind].ip))/2.5)*10.0^6
no_data=where(dat[ind].dip lt 0.0)
esed[no_data,4]=-99.0
esed[*,5]=(((alog(10))*(dat[ind].dzp)*(dat[ind].zp))/2.5)*10.0^6
no_data=where(dat[ind].dzp lt 0.0)
esed[no_data,5]=-99.0
esed[*,6]=(((alog(10))*(dat[ind].dIB427)*(dat[ind].IB427))/2.5)*10.0^6
no_data=where(dat[ind].dIB427 lt 0.0)
esed[no_data,6]=-99.0
esed[*,7]=(((alog(10))*(dat[ind].dIB464)*(dat[ind].IB464))/2.5)*10.0^6
no_data=where(dat[ind].dIB464 lt 0.0)
esed[no_data,7]=-99.0
esed[*,8]=(((alog(10))*(dat[ind].dIB484)*(dat[ind].IB484))/2.5)*10.0^6
no_data=where(dat[ind].dIB484 lt 0.0)
esed[no_data,8]=-99.0
esed[*,9]=(((alog(10))*(dat[ind].dIB505)*(dat[ind].IB505))/2.5)*10.0^6
no_data=where(dat[ind].dIB505 lt 0.0)
esed[no_data,9]=-99.0
esed[*,10]=(((alog(10))*(dat[ind].dIB527)*(dat[ind].IB527))/2.5)*10.0^6
no_data=where(dat[ind].dIB527 lt 0.0)
esed[no_data,10]=-99.0
esed[*,11]=(((alog(10))*(dat[ind].dIB574)*(dat[ind].IB574))/2.5)*10.0^6
no_data=where(dat[ind].dIB574 lt 0.0)
esed[no_data,11]=-99.0
esed[*,12]=(((alog(10))*(dat[ind].dIB624)*(dat[ind].IB624))/2.5)*10.0^6
no_data=where(dat[ind].dIB624 lt 0.0)
esed[no_data,12]=-99.0
esed[*,13]=(((alog(10))*(dat[ind].dIB679)*(dat[ind].IB679))/2.5)*10.0^6
no_data=where(dat[ind].dIB679 lt 0.0)
esed[no_data,13]=-99.0
esed[*,14]=(((alog(10))*(dat[ind].dIB709)*(dat[ind].IB709))/2.5)*10.0^6
no_data=where(dat[ind].dIB709 lt 0.0)
esed[no_data,14]=-99.0
esed[*,15]=(((alog(10))*(dat[ind].dIB738)*(dat[ind].IB738))/2.5)*10.0^6
no_data=where(dat[ind].dIB738 lt 0.0)
esed[no_data,15]=-99.0
esed[*,16]=(((alog(10))*(dat[ind].dIB767)*(dat[ind].IB767))/2.5)*10.0^6
no_data=where(dat[ind].dIB767 lt 0.0)
esed[no_data,16]=-99.0
esed[*,17]=(((alog(10))*(dat[ind].dIB827)*(dat[ind].IB827))/2.5)*10.0^6
no_data=where(dat[ind].dIB827 lt 0.0)
esed[no_data,17]=-99.0
esed[*,18]=(((alog(10))*(dat[ind].dNB711)*(dat[ind].NB711))/2.5)*10.0^6
no_data=where(dat[ind].dNB711 lt 0.0)
esed[no_data,18]=-99.0
esed[*,19]=(((alog(10))*(dat[ind].dNB816)*(dat[ind].NB816))/2.5)*10.0^6
no_data=where(dat[ind].dNB816 lt 0.0)
esed[no_data,19]=-99.0
esed[*,20]=(((alog(10))*(dat[ind].dks)*(dat[ind].ks))/2.5)*10.0^6
no_data=where(dat[ind].dks lt 0.0)
esed[no_data,20]=-99.0
esed[*,21]=(((alog(10))*(dat[ind].du_s)*(dat[ind].u_s))/2.5)*10.0^6
no_data=where(dat[ind].du_s lt 0.0)
esed[no_data,21]=-99.0
esed[*,22]=(((alog(10))*(dat[ind].dg_s)*(dat[ind].g_s))/2.5)*10.0^6
no_data=where(dat[ind].dg_s lt 0.0)
esed[no_data,22]=-99.0
esed[*,23]=(((alog(10))*(dat[ind].dr_s)*(dat[ind].r_s))/2.5)*10.0^6
no_data=where(dat[ind].dr_s lt 0.0)
esed[no_data,23]=-99.0
esed[*,24]=(((alog(10))*(dat[ind].di_s)*(dat[ind].i_s))/2.5)*10.0^6
no_data=where(dat[ind].di_s lt 0.0)
esed[no_data,24]=-99.0
esed[*,25]=(((alog(10))*(dat[ind].dz_s)*(dat[ind].z_s))/2.5)*10.0^6
no_data=where(dat[ind].dz_s lt 0.0)
esed[no_data,25]=-99.0
esed[*,26]=(((alog(10))*(dat[ind].dF814W)*(dat[ind].F814W))/2.5)*10.0^6
no_data=where(dat[ind].dF814W lt 0.0)
esed[no_data,26]=-99.0
esed[*,27]=(dat[ind].FLUX_ERR)*10.0^3
no_data=where(dat[ind].FLUX_ERR lt 0.0)
esed[no_data,27]=-99.0
esed[*,28]=(dat[ind].et250)*10.0^3
no_data=where(dat[ind].et250 lt 0.0)
esed[no_data,28]=-99.0
esed[*,29]=(dat[ind].et350)*10.0^3
no_data=where(dat[ind].et350 lt 0.0)
esed[no_data,29]=-99.0
esed[*,30]=(dat[ind].et500)*10.0^3
no_data=where(dat[ind].et500 lt 0.0)
esed[no_data,30]=-99.0

;check which elements of sed array are good (i.e. have data)
sed_good=where(sed gt 0.0 );and sed ne 99.00,ngood)
esed_good=where(esed gt 0.0,engood)
sed_bad=where(sed le 0.0,nbad)
esed_bad=where(esed lt 0.0,engood)


match,sed_good,esed_good,suba,subb, COUNT = ngood
sed_good=sed_good[suba]

;-------------------THIS BIT IS EXTREMELY COMPLICATED--------------------
;-----------------ADVISE AGAINST EDITING------------------------------
;;----(took me weeks to get working)---------------

;define data structure used by NMF_sparse
data=create_struct('nx',nbands*nz,'ny',nsrc,'rowstart',fltarr(nsrc),'nxrow',fltarr(nsrc),'val',fltarr(ngood),'x',fltarr(ngood))
data_ivar=create_struct('nx',nbands*nz,'ny',nsrc,'rowstart',fltarr(nsrc),'nxrow',fltarr(nsrc),'val',fltarr(ngood),'x',fltarr(ngood))

;build data
;for each source
cum_ind=0.0
;print, nsrc
for i=0,nsrc-1 do begin
xz=dat[ind[i]].redshift
;interpolate
iz=round(interpol(findgen(nz),tred,xz))
;print, iz,xz
;iz is redshift index for xz,
sed_ind=where(sed[i,*] gt 0.0,ngood)
esed_ind=where(esed[i,*] gt 0.0,engood)
match,sed_ind,esed_ind,suba,subb, COUNT = nsed_ind
sed_ind=sed_ind[suba]

data.val[cum_ind:nsed_ind+cum_ind-1]=sed[i,sed_ind]
data.x[cum_ind:nsed_ind+cum_ind-1]=sed_ind+(iz*nbands)
data.rowstart[i]=cum_ind
data.nxrow[i]=n_elements(sed_ind)

data_ivar.val[cum_ind:nsed_ind+cum_ind-1]=1.0/esed[i,sed_ind]^2
data_ivar.x[cum_ind:nsed_ind+cum_ind-1]=sed_ind+(iz*nbands)
data_ivar.rowstart[i]=cum_ind
data_ivar.nxrow[i]=n_elements(sed_ind)

cum_ind=cum_ind+n_elements(sed_ind)
ENDFOR
save, data,data_ivar,nbands,sed,esed,filename='/export/scratch/rh273/saves/data_sparse.sav'
END

pro run_NMF_fls_z_lt_1
;-------------code for running NMF_sparse---------
;requires that models have already been convolved with filters
;and the data has been put into correct format


;-----load data and model grid--------------
restore,'~/idl/NMF/Rachels_Code_Bits/data_sparse.sav'
restore,'~/idl/NMF/Rachels_Code_Bits/irgrid_test.sav'
;------------reform model grid into 2d array-------------
mmatrix=fltarr(nz*nbands,n_elements(IRGRID[*,0,0]))


for i=0,n_elements(IRGRID[*,0,0])-1 DO mmatrix[*,i]=reform(IRGRID[i,*,*],nz*nbands)

data_ivar.val[8080]=100.0
;--------------run NMF sparse-----------------------------
;no_temps=number of templates
nmf_sparse,data,data_ivar,3,mmatrix,3000,temp=temps,coeffs=coeffs
save, mmatrix,temps,coeffs,data,data_ivar,filename='~/idl/NMF/Rachels_Code_Bits/nmf_run_3000.sav'
END

pro NMF_sparse_plot
restore,'/export/scratch/rh273/saves/nmf_run_3000.sav'
restore,'/export/scratch/rh273/saves/mmodels.sav'
 AGN=mrdfits('/export/scratch/rh273/saves/AGNmod.fits', 1)
  SB=mrdfits('/export/scratch/rh273/saves/SBmodels.fits',1)
  restore,'/export/scratch/rh273/saves/cirrus_models.dat'
    templates=mmatrix#temps
  no_templates=1
  compo=fltarr(3,n_elements(templates[0,*]))
  for j=0,no_templates-1 DO BEGIN
  AGN_L_IR=0
  for i=0,n_elements(AGN)-1 DO BEGIN
  AGN_L_IR=AGN_L_IR+int_tabulated2(lambdaGRID,temps[i,j]*mmatrix[*,i])
  ENDFOR
  SB_L_IR=0
  for i=n_elements(AGN),n_elements(AGN)+n_elements(SB)-1 DO BEGIN
  SB_L_IR=SB_L_IR+int_tabulated2(lambdaGRID,temps[i,j]*mmatrix[*,i])
  ENDFOR
   Cir_L_IR=0
  for i=n_elements(AGN)+n_elements(SB),n_elements(mmatrix[0,*])-1 DO BEGIN
  Cir_L_IR=Cir_L_IR+int_tabulated2(lambdaGRID,temps[i,j]*mmatrix[*,i])
  ENDFOR
 L_IR=int_tabulated2(lambdaGRID,templates[*,j])
compo[*,j]=[AGN_L_IR/L_IR,SB_L_IR/L_IR,Cir_L_IR/L_IR];, (AGN_L_IR+SB_L_IR+Cir_L_IR)/L_IR
 ENDFOR
set_plot, 'ps'
device, filename='/export/scratch/rh273/saves/new_sparse_SED_templates.ps', /col
colours
!p.font=0
!p.thick=2
!x.thick=3
!y.thick=3
!p.charsize=1.2
!p.psym=0
colarray=[fsc_color('brown'),fsc_color('blue'),fsc_color('red'), fsc_color('green'), fsc_color('orange'), fsc_color('violet')]
plot, lambdaGRID, templates[*,0], /xlog, /ylog, yrange=[10.0^(-11),10.0^(-1)], xtitle=textoidl('\lambda (micron)') , ytitle=('Flux scaled');, xTICKFORMAT='(E10.4)'
xpos=[10,60,1000]
comps=['AGN %','Starburst %','Cirrus %']
for j=0,2 DO XYOUTS,xpos[j],10.0^(-1.0*(6)),comps[j]
xpos=[10,100,1000]

for i=0,n_elements(templates[0,*])-1 DO BEGIN
oplot, lambdaGRID, templates[*,i], color=colarray[i]
for j=0,2 DO BEGIN
comps=string(compo[j,i],format='(F5.2)')
XYOUTS,xpos[j],10.0^(-1.0*((i+1)*0.7+6)),comps,color=colarray[i]
ENDFOR

ENDFOR

CLOSEPS
END





pro NMF_template_composition
restore,'/export/scratch/rh273/saves/nmf_run.sav'
restore,'/export/scratch/rh273/saves/mmodels.sav'

  AGN=mrdfits('/export/scratch/rh273/saves/AGNmod.fits', 1)
  SB=mrdfits('/export/scratch/rh273/saves/SBmodels.fits',1)
  restore,'/export/scratch/rh273/saves/cirrus_models.dat'
  
  
  templates=mmatrix#temps
  for j=0,5 DO BEGIN
  AGN_L_IR=0
  for i=0,n_elements(AGN)-1 DO BEGIN
  AGN_L_IR=AGN_L_IR+int_tabulated2(lambdaGRID,temps[i,j]*mmatrix[*,i])
  ENDFOR
  SB_L_IR=0
  for i=n_elements(AGN),n_elements(AGN)+n_elements(SB)-1 DO BEGIN
  SB_L_IR=SB_L_IR+int_tabulated2(lambdaGRID,temps[i,j]*mmatrix[*,i])
  ENDFOR
   Cir_L_IR=0
  for i=n_elements(AGN)+n_elements(SB),n_elements(mmatrix[0,*])-1 DO BEGIN
  Cir_L_IR=Cir_L_IR+int_tabulated2(lambdaGRID,temps[i,j]*mmatrix[*,i])
  ENDFOR
 L_IR=int_tabulated2(lambdaGRID,templates[*,j])
 print, AGN_L_IR/L_IR,SB_L_IR/L_IR,Cir_L_IR/L_IR, (AGN_L_IR+SB_L_IR+Cir_L_IR)/L_IR
 ENDFOR
  set_plot, 'ps'
device, filename='/export/scratch/rh273/saves/SED_templates_comp.ps', /col
colours
!p.font=1
!p.thick=2
!x.thick=3
!y.thick=3
!p.charsize=1.2
!p.psym=0
plot, temps[*,0], /ylog
colarray=[fsc_color('blue'),fsc_color('red'), fsc_color('green'), fsc_color('orange'), fsc_color('violet')]

for i=1,5 DO oplot, temps[*,i], color=colarray[i-1]
oplot,[n_elements(AGN)-1,n_elements(AGN)-1],[min(temps),1], linestyle=1
oplot,[n_elements(AGN)+n_elements(SB)-1,n_elements(AGN)+n_elements(SB)-1],[min(temps),1], linestyle=2

CLOSEPS
END

pro templates_fit_FLS_z_lt_1
;convolve NMF derived templates with filters

 infilt=[137, 164, 165, 166, 167, 168, 121]
dir='/export/scratch/rh273/saves/'
filt='filter_seb'

tred=findgen(200)/200.0
nz=200
nif=n_elements(infilt)  ;number of infrared filters
fiwav=fltarr(nif,900) ;lambda array for all infra-red filters
fitran=fltarr(nif,900)  ;transmission array for all infrared filters


;load NMF-ir temps
restore,'/export/scratch/rh273/saves/nmf_run.sav'
restore,'/export/scratch/rh273/saves/mmodels.sav'
;uses mmatrix from mmodels
  templates=mmatrix#temps


nit=n_elements(templates[0,*]); number of templates
iwav=fltarr(nit,n_elements(lambdaGRID))
for i=0,nit-1 DO iwav[i,*]=lambdaGRID*10^4 ; lambda grid for each model, in angstroms
itflx=transpose(templates)
for i=0,nit-1 DO itflx[i,*]=(itflx[i,*]*3*10.0^18)/iwav[i,*]^2 ; convert to SED flux density per unit wavelength (as required by filt_mag)

;load filters
irfeff=fltarr(nif) ; effective lambda for all infrared filters
irfsurf=fltarr(nif) ; surface function for all infrared filters
for k=0,n_elements(infilt)-1 do begin
;uses Seb's filt_read, loops for one res file, over filter ids
filt_read,filt,dir=dir,infilt[k],np,name,f1wav,f1tran
fiwav[k,0:np-1]=f1wav[0:np-1] 
fitran[k,0:np-1]=f1tran[0:np-1]/max(f1tran); normalised transmission
filt_properties,f1wav,f1tran,eff_lam=teff,surf=tsurf ;calculate effective lambda and surface function
irfeff[k]=teff
irfsurf[k]=tsurf
endfor

itmp=intarr(nit) ; number of infra red templates
irgrid=fltarr(nit,nif,nz) ; array number of templates by number of filters, by number of redshifts
 
 for i=0,nit-1 do begin
 for ii=0,nif-1 do begin
 for iii=0,nz-1 do begin
tz=tred[iii]
red=tred[iii]
tpts=where(iwav[i,*] ne 0,ntpts) ; where lambda template /= 0
;creates arrays that exclude values for when lambda =0
tiwav=fltarr(ntpts)
tiflx=fltarr(ntpts)
fpts=where(fiwav[ii,*] ne 0,nfpts)
tfwav=fltarr(nfpts)
tftran=fltarr(nfpts)
;fills arrays for values for when lambda /=0
tiwav[0:ntpts-1]=iwav[i,tpts]
tiflx[0:ntpts-1]=[itflx[i,tpts]]
tfwav[0:nfpts-1]=fiwav[ii,fpts]/(1.+tz) ;move filter to actual emission wavelength
tftran[0:nfpts-1]=fitran[ii,fpts]


filt_mag,tiwav,tiflx,tfwav,tftran,tmp;get magnitude for template i through filter ii

irgrid[i,ii,iii]=10^((tmp)*(-0.4))/irfsurf[ii]*(irfeff[ii]/(1.+red))^2/2.99792e18;convert flux(freq)
 
;tred[k]=tz
itmp[i]=i
;print,tred(k),cold_kcor(k),sb9_kcor(k)
endfor
endfor
ENDFOR
save, irgrid,tred,nz, filename='/export/scratch/rh273/saves/templates_grid.sav'
END
