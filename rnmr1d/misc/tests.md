## Usage

    $ docker run -it --rm rnmr1d -h

```
Loading required package: methods
Loading required package: foreach
Loading required package: iterators
Loading required package: parallel
Rnmr1D.

Usage:
  Rnmr1D [options]

Options:
  -h, --help           Show this screen.
  -d, --debug          Show more information
  -z, --zip <file>     the full path name of the ZIP file (raw.zip)
  -s, --samples <file> the full path name of the Sample file (tabular format)
  -p, --proccmd <file> the full path name of the Macro-commands file for processing (text format)
  -b, --bucfile <file> the full path name of the file of bucket's zones (tabular format)
  -c, --cpu <n>        the number of cores
  -o, --outdir <path>  the full path name of the directory to output the resulting files [default: /data]
  -l, --logfile <file> the full path name of the LOG file [default: stderr] 
---
```

## Run the provided example script './examples/scripts/script_test.sh'

    $ sh ./examples/scripts/script_test.sh 

```
Loading required package: methods
Loading required package: foreach
Loading required package: iterators
Loading required package: parallel
[1] 1
Rnmr1D:  --- READING and CONVERTING ---
[1/23]: BPA_c21_aq_133-BPA250ng
-----
Read the 1R ...

[2/23]: BPA_c21_aq_135-BPA250ng
-----
Read the 1R ...

[3/23]: BPA_c21_aq_136-BPA25ng
-----
Read the 1R ...

[4/23]: BPA_c21_aq_137-BPA25ng
-----
Read the 1R ...

[5/23]: BPA_c21_aq_138-1-BPA250ng
-----
Read the 1R ...

[6/23]: BPA_c21_aq_138-2-BPA250ng
-----
Read the 1R ...

[7/23]: BPA_c21_aq_138-3-BPA250ng
-----
Read the 1R ...

[8/23]: BPA_c21_aq_141-BPA250ng
-----
Read the 1R ...

[9/23]: BPA_c21_aq_149-BPA25ng
-----
Read the 1R ...

[10/23]: BPA_c21_aq_152-BPA25ng
-----
Read the 1R ...

[11/23]: BPA_c21_aq_154-BPA250ng
-----
Read the 1R ...

[12/23]: BPA_c21_aq_155-BPA250ng
-----
Read the 1R ...

[13/23]: BPA_c21_aq_157-BPA250ng
-----
Read the 1R ...

[14/23]: BPA_c21_aq_162-BPA25ng
-----
Read the 1R ...

[15/23]: BPA_c21_aq_163-BPA25ng
-----
Read the 1R ...

[16/23]: BPA_c21_aq_164-BPA25ng
-----
Read the 1R ...

[17/23]: BPA_c21_aq_165-BPA250ng
-----
Read the 1R ...

[18/23]: BPA_c21_aq_166-BPA25ng
-----
Read the 1R ...

[19/23]: BPA_c21_aq_172-BPA250ng
-----
Read the 1R ...

[20/23]: BPA_c21_aq_173-BPA250ng
-----
Read the 1R ...

[21/23]: BPA_c21_aq_177-BPA25ng
-----
Read the 1R ...

[22/23]: BPA_c21_aq_179-BPA250ng
-----
Read the 1R ...

[23/23]: BPA_c21_aq_181-BPA25ng
-----
Read the 1R ...

Rnmr1D:  Generate the final matrix of spectra...
Rnmr1D:  Write the spec.pack file ...
Rnmr1D:  Write the list_pars.csv file ...
Rnmr1D: ------------------------------------
Rnmr1D: Process the Macro-commands file
Rnmr1D: ------------------------------------
Rnmr1D: 
Rnmr1D:  Calibration: PPM REF =0, Zone Ref = (-0.013,0.013)
Rnmr1D:  Baseline Correction: PPM Range = ( -0.4995660520046 , 11.0000011357711 )
Rnmr1D:  Window Size = 256
Rnmr1D:  Baseline Correction: PPM Range = ( 6.781 , 7.406 )
Rnmr1D:  Window Size = 26
Rnmr1D:  Alignment of the selected zones ...
Rnmr1D:  Zone 1 = ( 7.758 , 7.809 )
Rnmr1D:  Zone 2 = ( 7.048 , 7.083 )
Rnmr1D:  Zone 3 = ( 6.969 , 7.045 )
Rnmr1D:  Zeroing the selected PPM ranges ...
Rnmr1D:  Zone 1 = ( 6.783 , 6.879 )
Rnmr1D:  Zone 2 = ( 7.081 , 7.174 )
Rnmr1D:  Zone 3 = ( 6.928 , 6.989 )
Rnmr1D:  Alignment of the selected zones ...
Rnmr1D:  Zone 1 = ( 8.103 , 8.157 )
Rnmr1D:  Baseline Correction: PPM Range = ( 8.158 , 8.42 )
Rnmr1D:  Window Size = 26
Rnmr1D:  Baseline Correction: PPM Range = ( 5.882 , 6.178 )
Rnmr1D:  Window Size = 26
Rnmr1D:  Zeroing the selected PPM ranges ...
Rnmr1D:  Zone 1 = ( 4.536 , 5.094 )
Rnmr1D:  Baseline Correction: PPM Range = ( 2.62 , 4.561 )
Rnmr1D:  Window Size = 43
Rnmr1D:  Baseline Correction: PPM Range = ( 1.981 , 2.61 )
Rnmr1D:  Window Size = 43
Rnmr1D:  Baseline Correction: PPM Range = ( 1.573 , 1.985 )
Rnmr1D:  Window Size = 43
Rnmr1D:  Baseline Correction: PPM Range = ( 1.296 , 1.617 )
Rnmr1D:  Window Size = 26
Rnmr1D:  Baseline Correction: PPM Range = ( 1.128 , 1.297 )
Rnmr1D:  Window Size = 26
Rnmr1D:  Baseline Correction: PPM Range = ( 0.824 , 1.092 )
Rnmr1D:  Window Size = 26
Rnmr1D:  Zeroing the selected PPM ranges ...
Rnmr1D:  Zone 1 = ( 0.814 , 0.897 )
Rnmr1D:  Alignment of the selected zones ...
Rnmr1D:  Zone 1 = ( 8.6 , 8.64 )
Rnmr1D:  Zone 2 = ( 8.573 , 8.6 )
Rnmr1D:  Zone 3 = ( 8.327 , 8.362 )
Rnmr1D:  Baseline Correction: PPM Range = ( 5.584 , 5.663 )
Rnmr1D:  Window Size = 26
Rnmr1D:  Alignment of the selected zones ...
Rnmr1D:  Zone 1 = ( 3.401 , 3.455 )
Rnmr1D:  Zone 2 = ( 3.243 , 3.287 )
Rnmr1D:  Baseline Correction: PPM Range = ( 3.753 , 3.807 )
Rnmr1D:  Window Size = 26
Rnmr1D:  Alignment of the selected zones ...
Rnmr1D:  Zone 1 = ( 3.779 , 3.797 )
Rnmr1D:  Alignment of the selected zones ...
Rnmr1D:  Zone 1 = ( 3.753 , 3.779 )
Rnmr1D:  Alignment of the selected zones ...
Rnmr1D:  Zone 1 = ( 3.71 , 3.753 )
Rnmr1D:  Write the spec.pack file ...
Rnmr1D: 
Rnmr1D: ------------------------------------
Rnmr1D: Process the file of buckets
Rnmr1D: ------------------------------------
Rnmr1D: 
Rnmr1D: NB Buckets = 211
Rnmr1D: 
Rnmr1D: 
   user  system elapsed 
  3.013   1.063   2.741 
```
