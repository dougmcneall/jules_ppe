# Constraining the historic carbon cycle in JULES-ES1.0
D.J. McNeall, E. Robertson & A. Wiltshire

Code and data for the generation and analysis of a perturbed parameter ensemble of land surface model JULES-ES-1.0
D.McNeall 10th October 2022  


## Ensemble design and generation
Ensemble generation proceeded in two waves, with increasingly organised code:

Wave00 first 300 members generated by write_design_u-ao732.R, which calls write_jules_design2.R and writes JULES 
configuration files to folder /conf_files_u-ao732.

Further 200 members generated by write_design_u-ao732a.R and standard member generated by write_standard_u-ao732.R
which write configuration files to folder /conf_files_u-ao732a.

These configuration files were then passed to a Rose suite for running the model. Global summary data can be found in the /data folder.

Second wave generated by design-JULES-ES-1p0-wave01.Rmd, which writes configuration files to folder 
/conf_files_augment_JULES-ES-1p0.


The constraint part of the analysis is coded in constrain-JULES-ES-1p0.Rmd, the sensitivity analysis in 
sensitivity-JULES-ES-1p0.Rmd and emulator testing and validation is in test-emulators-JULES-ES-1p0.Rmd.







