# Overview

This folder contains supporting material for MD simulation results presented in

### Ovchinnikov, Louveau, Barton, Karplus and Chakraborty

Matlab data files and scripts for generating the figures are described below.

# File list and descriptions

`3bnc60at-rmsf-cg-hc.mat`  - RMSF datafile for intermediate 3BNC60 HC  
`3bnc60at-rmsf-cg-lc.mat`  - RMSF datafile for intermediate 3BNC60 LC  
`3bnc60glt-rmsf-cg-hc.mat` - RMSF datafile for germline 3BNC60 HC  
`3bnc60glt-rmsf-cg-lc.mat` - RMSF datafile for germline 3BNC60 LC  
`3bnc60t-rmsf-cg-hc.mat`   - RMSF datafile for mature 3BNC60 HC  
`3bnc60t-rmsf-cg-lc.mat`   - RMSF datafile for mature 3BNC60 LC  

`3h109lt-rmsf-cg-hc.mat` - RMSF datafile for intermediate PGT121 HC  
`3h109lt-rmsf-cg-lc.mat` - RMSF datafile for intermediate PGT121 LC  
`gl121t-rmsf-cg-hc.mat`  - RMSF datafile for germline PGT121 HC  
`gl121t-rmsf-cg-lc.mat`  - RMSF datafile for germline PGT121 LC  
`pgt121t-rmsf-cg-hc.mat` - RMSF datafile for mature PGT121 HC  
`pgt121t-rmsf-cg-lc.mat` - RMSF datafile for mature PGT121 LC  

`ch103-i3.2t-rmsf-cg-hc.mat` - RMSF datafile for intermediate CH103 HC  
`ch103-i3.2t-rmsf-cg-lc.mat` - RMSF datafile for intermediate CH103 LC  
`ch103t-rmsf-cg-hc.mat`      - RMSF datafile for mature CH103 HC  
`ch103t-rmsf-cg-lc.mat`      - RMSF datafile for mature CH103 LC  
`ch103ucat-rmsf-cg-hc.mat`   - RMSF datafile for germline CH103 HC  
`ch103ucat-rmsf-cg-lc.mat`   - RMSF datafile for germline CH103 LC  

`malign-lc.mat` - sequence alignment file for heavy chains  
`malign.mat` - sequence alignment file for heavy chains  

`entropy-ca1.mat` - datafile that contains quasiharmonic entropies computed from MD (1st trajectory only)  
`entropy-ca2.mat` - datafile that contains quasiharmonic entropies computed from MD (2nd trajectory only)  
`entropy-ca3.mat` - datafile that contains quasiharmonic entropies computed from MD (3rd trajectory only)  
`entropy-ca4.mat` - datafile that contains quasiharmonic entropies computed from MD (4th trajectory only)  
`entropy-ca5.mat` - datafile that contains quasiharmonic entropies computed from MD (5th trajectory only)  

`colours.m` - aux matlab file to define plot colors  
`showdom.m` - auxiliary matlab file for plotting  
`smooth2.m` - aux matlab file for plotting  

`mkrms.m` - matlab main file to generate RMSF figures  

### Usage

`*.mat` files are in standard matlab data format; they can be loaded into matlab using `load(<filename>)`  
`*.m` files are matlab scripts; they can be run using `run <filename>`  

for example:  
`run mkrms` will load the RMSF data computed from MD and regenerate the RMSF plots in the paper (Fig. 1)  
`run mkentropy` will load the quasiharmonic entropy data computed from MD and regenerate the entropy plots in the paper (Fig. 3)  

Matlab version R2010b was used to prepare the figures
