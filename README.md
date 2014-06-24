#README
------------------
##Contents
- filter_master.*: filter files created with code in ./filters 	
- ./filters: Directory containing Seb's filter codes and data. Code required by templates_fit_fls_z_lt_1.pro

- mmodels.sav : Michael and Andreas's SB, Cirrus and AGN models all in one matrix
- templates_fit_fls_z_lt_1.pro: idl code. Contains numerous programs. 
  1. First runs models through desired filters,
  2. Second puts data into required sparse structure
  3. Third runs the NMF_sparse routine (requires idlutils lib to be installed)