PRO FILT_read, file_stem, dir=dir, filter_id, npoint, name, wave, trans, $
 silent=silent
;, $
;               eff_lam=eff_lam, surf=surf, band=band

;+
; NAME:
;	FILT_read
;
; PURPOSE:
;
;	reads a specfic filter transmission from a FILTER.RES file (used in 
;       e.g. HYPERZ)
;
;
; CATEGORY:
;	analysis
;
; CALLING SEQUENCE:
;
;	filt_read, file_stem, dir=dir, filter_id, npoint, name, wave, trans,  silent=silent
;
;
; INPUTS:
;
;      file_stem: name of .res file (without the ".res")
;      filter_id: number of filter to be read
;      	
;
; OPTIONAL INPUTS:
;	
;	
; KEYWORD PARAMETERS:
;       dir:  directory containing filter files (default ".")
;	silent: if set Runs without output
;
; OUTPUTS:
;	npoint:  number of data points in filter
;       name:    name of filter
;       wave:    Wavelength in AA
;       trans:   Transmision
;
; OPTIONAL OUTPUTS:
;	
;
; COMMON BLOCKS:
;	
;
; SIDE EFFECTS:
;	
;
; RESTRICTIONS:
;
;  N.B. the correct .log file must have been generated using
;  FILT_MAKE_DAT	
;
; PROCEDURE:
;	
;
; EXAMPLE:
;	
;
; MODIFICATION HISTORY:
; 	Written by:	Seb Oliver (Sussex) 19th Mar 2000
;	July, 1996	
;-

   IF n_params() NE 6 THEN message, ' CALLING SEQUENCE: filt_read, file_stem, dir=dir, filter_id, npoint, name, wave, trans,  silent=silent'

   IF NOT keyword_set(silent) THEN silent = 0
   IF NOT keyword_set(dir) THEN dir = '.' 
   file = findfile( dir+'/'+strcompress(file_stem, /rem)+'.log')
   dat_file = findfile( dir+'/'+strcompress(file_stem, /rem)+'.dat')
   IF file[0] EQ '' THEN BEGIN
      message, 'Generating summary file', /inf
      filt_make_dat, file_stem, dir=dir
   ENDIF

; no need for a skip line as header lines generate an error and are ignored   
   readfmt, dir+'/'+strcompress(file_stem, /rem)+'.log', '(i6,x,a63,i4)',  id, name, npoints,  silent=silent

; various checks
;
; check ID numbers are sequential
   check_ids = id-indgen(n_elements(id))-1
   tmp = where(check_ids ne 0, nbad)
   IF nbad GT 0 THEN BEGIN 
      message, 'ERROR in some filter ids', /inf
      print, name[tmp]
      print, id[tmp]
      message, 'ABORTING'
   ENDIF

; check number of points adds up to number of lines in .res file

   nlines = numlines( dir+'/'+strcompress(file_stem, /rem)+'.res')
   check_numlines = round(total(npoints))+n_elements(id)
   IF  check_numlines NE nlines THEN message, 'Number of lines do not match expected'


; selecting relevent filter

   select = where(id EQ filter_id, n_id)
   IF n_id EQ 0 OR n_elements(select) GT 1 THEN message, 'Requested Filter not present'
   select = select[0]

; print selected details
   IF silent eq 0 THEN BEGIN 
      print, 'Selecting:'
      print, 'ID:', id[select]
      print, 'Name:', name[select]
      print, 'Npoints', npoints[select]
   ENDIF


; REading res file
   get_lun, unit
   openr, unit, dir+'/'+strcompress(file_stem, /rem)+'.res'
   
;skiping some lines
   IF select GT 0 THEN  skip = round(total(npoints(0:(select-1))+1)) ELSE skip = 0
;   stop
   a = 'a'
   FOR i=0L, skip-1L DO readf, unit, a

   FILT_read_next, unit, npoint, name, wave, trans

   IF unit GT 0 THEN free_lun, unit

   IF silent EQ 0 THEN BEGIN 
      print, 'Read:'
      print, 'Name:', name
      print, 'Npoints', npoint
   ENDIF


   
END
