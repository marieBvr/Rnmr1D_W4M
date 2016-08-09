# Rnmr1D

## Installation

Clone this repository, then `cd` to your clone path.

    $ git clone https://bitbucket.org/nmrprocflow/rnmr1d
    $ cd rnmr1d

### Standard R package

Requirements:

    R packages: Rcpp (>= 0.11.3), docopt (>= 0.4.2), stringr (>= 1.0.0), doParallel (>= 1.0.8)

Install this package like any R standard package:

    $ R CMD INSTALL .

Copy the 'exec/Rnmr1D' file to '/usr/local/bin', then test it:

    $ cp ./exec/Rnmr1D /usr/local/bin
    $ which Rnmr1D
    $ Rnmr1D -h

Run the provided example script:

```
    $ sh ./examples/scripts/test.sh 
```

### Docker

Create a docker container to use this library via:

    $ docker build -t rnmr1d .

Usage:

    $ docker run -it --rm rnmr1d -h

```
Rnmr1D - Command Line Interface (CLI) of the NMR spectra processing module (R package 'Rnmr1D')

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
  -l, --logfile <file> the full path name of the LOG file [default: stderr]
```


Run the provided example script:

```
    $ sh ./examples/scripts/docker_test.sh 
```

