
pro filter_make_master, files, dir, file_output
;+
; NAME:
;
;
;
; PURPOSE:
;Make a master res file from .res files
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;filter_make_master,files,dir
;
;
; INPUTS:
;files = string array of filter files to be used in master file
;dir = directory of the filter files (all files are asssumed to be in
;same directory)
;file_output = filename for output .res, .log .sum etc files
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;creates . res, .sum , .log files
;
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;created by Peter Hurley 31/07/2012
;-
command='cat '
for i=0,n_elements(files)-1 DO command=command+files[i]+'.res '
command=command+'>'+file_output+'.res'
spawn, command
filt_make_dat,file_output,dir=dir
END
