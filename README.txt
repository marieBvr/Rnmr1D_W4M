## ****** Determine Rnmr1D environnemnt : ****** ## 
# Version 2016 
# M Lefebvre / D Jacob
# Developed within MetaboHUB, The French Infrastructure in Metabolomics and Fluxomics (www.metabohub.fr/en)

## --- PERL compilator / libraries : --- ## 
NA 
-- 

## --- R bin and Packages : --- ## 
R version 3.1.1 (2014-07-10) -- "Sock it to Me"
Copyright (C) 2014 The R Foundation for Statistical Computing
Platform: x86_64-redhat-linux-gnu (64-bit)

 # Install the R libraries, required for Rnmr1D #
 !WARNING! : You should respect installation order ! 
   # You can find binaries in the tool folder /static/binaries #
   R CMD INSTALL ./iterators_1.0.8.tar.gz
   R CMD INSTALL ./foreach_1.4.3.tar.gz
   R CMD INSTALL ./doParallel_1.0.10.tar.gz
   R CMD INSTALL ./stringi_1.0-1.tar.gz
   R CMD INSTALL ./magrittr_1.5.tar.gz
   R CMD INSTALL ./stringr_1.0.0.tar.gz
   R CMD INSTALL ./docopt_0.4.3.3.tar.gz
   R CMD INSTALL ./Rcpp_0.12.3.tar.gz

   # Install Rnmr1D #
    # In the tool folder #
    R CMD INSTALL .
    cp ./exec/Rnmr1D /usr/local/bin
    Rnmr1D -h
   # The lastest command should display : #
   "Rnmr1D - Command Line Interface (CLI) of the NMR spectra processing module (R package 'Rnmr1D')

   Usage:
    Rnmr1D [options]

    Options:
      -h, --help           Show this screen.
      -d, --debug          Show more information
      -z, --zip <file>     the full path name of the ZIP file (raw.zip)
      -s, --samples <file> the full path name of the Sample file (tabular format)
      -p, --proccmd <file> the full path name of the Macro-commands file for processing (text format)
      -b, --bucfile <file> the full path name of the file of bucket's zones (tabular format)
      -c, --cpu <n>        the number of cores [default: 4]
      -o, --outdir <path>  the full path name of the directory to output the resulting files
      -l, --logfile <file> the full path name of the LOG file [default: stderr]"
--

## --- Binary dependencies --- ## 
NA 
-- 

## --- Config --- ## 
NA 
-- 

## --- XML HELP PART --- ## 
# -- images: 
Rnmr1D_NMRProcFlow_logo.png 
Rnmr1D-workflowPosition.png

# -- datatype :
# You need to install the datatypes "no_unzip.zip" for the Rnmr1D package.
Add the datatype to the $galaxybase/datatypes_conf.xml:
 <datatype extension="no_unzip.zip" type="galaxy.datatypes.no_unzip_datatypes:NoUnzip" display_in_upload="true" />	
Then copy the no_unzip_datatypes.py to your $galaxybase/lib/galaxy/datatypes/
-- 

## --- DATASETS --- ## 
No data set ! Waiting for galaxy pages 
-- 

