#------------------------------------------------
# Rnmr1D package: Build 1r spectrum from FID file (Bruker/Varian/nmrML)
# Project: NMRProcFlow - MetaboHUB
# (C) 2015 - D. JACOB - IMRA UMR1332 BAP
#------------------------------------------------

require(XML)
require(base64enc)

Spec1r.doProc <- function(Input, ...) UseMethod("Spec1r.doProc")

Spec1r.doProc.default <- function(Input, ...)
{
   .CALL(Input,...)
}

### Processing parameters
Spec1r.Procpar.default <- list (

### General Parameters
    DEBUG=TRUE,							# Debug 
    LOGFILE=stderr(),					# Messages output file
    VENDOR="bruker",					# Instrumental origin of the raw data, (bruker, varian)
    READ_RAW_ONLY=FALSE,				# Read Raw file only; do not carry out processing; if raw file is depending on INPUT_SIGNAL
    INPUT_SIGNAL="fid",					# What type of input signal: 'fid' or '1r'
    PDATA_DIR='pdata/1',				# subdirectory containing the 1r file (bruker's format only)

### Calibration
    TSP=TRUE,							# indicate if TSP signal present
    TSP_PPM_RANGE=1.5,					# indicate the ppm range where the TSP signal could be present if required
    PLASMA=FALSE,						# indicate if Sample type is plasma
    PPM_CALIBRATION=FALSE,				# indicate if PPM calibration is applied
    ZONE_REF=NULL,						# ppm range for calibration
    PPM_REF=0,							# ppm value of the center of resonance 

### PRE-PROCESSING
    LB= 0.3,							# Line Broadening parameter
    ZEROFILLING=TRUE,					# Zero Filling
    LINEBROADENING=TRUE,				# Line Broading

### Baseline correction parameters
    BASELINECOR=FALSE,					# indicate if Baseline correction is applied
    WINDOWSIZE=50,						# window size for the smoothing (spectral points)
    NEIGHSIZE=30,						# number of needed consecutive points to be eligible as baseline points (spectral points)
    NOISELEVEL=1.5,						# Factor of noise level

### Denoising
    DENOISING=FALSE

)

#--------------------------------
# Initialize Parameter Lists by the default ones
#--------------------------------
Spec1r.Procpar    <- Spec1r.Procpar.default

#--------------------------------
# verbose function
#--------------------------------
.v <- function(..., logfile=Spec1r.Procpar$LOGFILE) cat(sprintf(...), sep='', file=logfile, append=TRUE)

#--------------------------------
# READ FID && Acquisition Parameters
#--------------------------------

#### Retrieve a parameter value from Varian acquisition parameters
#--  internal routine
#--  ACQ: list of acquisition parameters of the acqus file
#--  paramStr: name of the parameter
.bruker.get_param <- function (ACQ,paramStr,type="numeric", arrayType=FALSE)
{
   regexpStr <- paste("^...",paramStr,"=",sep="")
   if (!arrayType) {
       acqval <- gsub("^[^=]+= ","" ,ACQ[which(simplify2array(regexpr(regexpStr,ACQ))>0)])
   } else {
       acqval <-simplify2array(strsplit(ACQ[which(simplify2array(regexpr(regexpStr,ACQ))>0)+1]," "))[2]
   }
   if (type=="numeric")
       acqval <-as.numeric(acqval)
   if (type=="string")
       acqval <-gsub("<","",gsub(">","",acqval))
   acqval
}

#### Retrieve a parameter value from Bruker acquisition parameters
#--  internal routine
#--  ACQ: list of acquisition parameters of the procpar file
#--  paramStr: name of the parameter
.varian.get_param <- function (ACQ,paramStr,type="numeric")
{
   regexpStr <- paste("^",paramStr," ",sep="")
   acqval <- gsub("1 ","",ACQ[which(simplify2array(regexpr(regexpStr,ACQ))>0)+1])
   if (type=="numeric")
       acqval <-as.numeric(acqval)
   if (type=="string")
       acqval <-gsub("\"","",acqval)
   acqval
}

#### Read FID data and the main parameters needed to generate the real spectrum
#--  external routine
#-- DIR: bruker directory containing the FID
.read.FID.bruker <- function(DIR)
{
   cur_dir <- getwd()
   setwd(DIR)

   # FID filename
   FIDFILE <- "fid"
   if (!file.exists(FIDFILE)) 
       stop("File ", FIDFILE, " does not exist\n")
   # Acquisition parameters filename
   ACQFILE <- "acqus"
   if (!file.exists(ACQFILE)) 
       stop("Acquisition parameters File (acqus) does not exist\n")

   ACQ     <- readLines(ACQFILE)
   PROBE   <- .bruker.get_param(ACQ,"PROBHD",type="string")
   SOLVENT <- .bruker.get_param(ACQ,"SOLVENT",type="string")
   PULSE   <- .bruker.get_param(ACQ,"PULPROG",type="string")
   TEMP    <- .bruker.get_param(ACQ,"TE",type="string")
   RELAXDELAY  <- .bruker.get_param(ACQ,"D",type="string",  arrayType=TRUE)
   PULSEWIDTH <- .bruker.get_param(ACQ,"P",type="string",  arrayType=TRUE)
   SPINNINGRATE <- .bruker.get_param(ACQ,"MASR",type="string")
   TD      <- .bruker.get_param(ACQ,"TD")
   SW      <- .bruker.get_param(ACQ,"SW")
   SWH     <- .bruker.get_param(ACQ,"SW_h")
   SFO1    <- .bruker.get_param(ACQ,"SFO1")
   GRPDLY  <- .bruker.get_param(ACQ,"GRPDLY")
   DTYPA   <- .bruker.get_param(ACQ,"DTYPA")
   BYTORDA <- .bruker.get_param(ACQ,"BYTORDA")
   ENDIAN = ifelse( BYTORDA==0, "little", "big")
   INSTRUMENT <- "Bruker"
   SOFTWARE <- gsub("##TITLE= Parameter file, ","", gsub("\t\tVersion"," v.",ACQ[1]))
   ORIGIN   <- gsub("^[^=]+= ","", ACQ[which(simplify2array(regexpr("^..ORIGIN=",ACQ))>0)])
   ORIGPATH <- gsub("^.. ", "", ACQ[which(simplify2array(regexpr("acqus$",ACQ))>0)])

   SIZE = ifelse( DTYPA==0, 4, 8)

   to.read = file(FIDFILE,"rb")
   signal<-readBin(to.read, what="integer",size=SIZE, n=TD, signed = TRUE, endian = ENDIAN)
   close(to.read)
   setwd(cur_dir)

   td <- length(signal)
   rawR <- signal[seq(from = 1, to = td, by = 2)]
   rawI <- signal[seq(from = 2, to = td, by = 2)]

   mediar<-mean(as.integer(rawR[c((3*length(rawR)/4):length(rawR))]),na.rm = TRUE)
   mediai<--mean(as.integer(rawI[c((3*length(rawR)/4):length(rawR))]),na.rm = TRUE)
   rawR<-rawR-mediar
   rawI<-rawI-mediai
   fid <- rawR+1i*rawI

   GRPDLY <- NULL
   GRPDLY <- tryCatch ( { if (!is.na(as.numeric(GRPDLY)) && as.numeric(GRPDLY)>0) round(GRPDLY); }, warning = function(w){ GRPDLY <- NULL;}, error=function(e){ GRPDLY <- NULL;} )
   if (is.null(GRPDLY)) {
      ### Estimation of point number of the dead time for the FID (Group Delay)
      nd0 <-which.max(abs(Re(fid)))
      nd0p <- which.max(abs(Re(fid)[1:(nd0-1)]))
      if ((abs(Re(fid)[nd0p])/abs(Re(fid)[nd0]))>0.5) nd0 <- nd0p;
      while(sign(Re(fid)[nd0-1])==sign(Re(fid)[nd0])) nd0 <- nd0-1;
      nd0 <- nd0-1
      GRPDLY <- nd0
   }
   proc <- list( INSTRUMENT=INSTRUMENT, SOFTWARE=SOFTWARE, ORIGIN=ORIGIN, ORIGPATH=ORIGPATH, 
                PROBE=PROBE, PULSE=PULSE, SOLVENT=SOLVENT, TEMP=TEMP, 
                RELAXDELAY=RELAXDELAY, SPINNINGRATE=SPINNINGRATE, PULSEWIDTH=PULSEWIDTH,
                TD=TD, SW=SW, SWH=SWH, SFO1=SFO1, GRPDLY=GRPDLY, phc0=0, phc1=0 )
   spec <- list( path=DIR, proc=proc, fid=fid )

   spec
}

.read.FID.varian <- function(DIR)
{
   cur_dir <- getwd()
   setwd(DIR)

   # FID filename
   FIDFILE <- "fid"
   if (!file.exists(FIDFILE)) 
       stop("File ", FIDFILE, " does not exist\n")
   # Acquisition parameters filename
   ACQFILE <- "procpar"
   if (!file.exists(ACQFILE)) 
       stop("Acquisition parameters File (procpar) does not exist\n")

   ACQ     <- readLines(ACQFILE)
   SWH <-  .varian.get_param(ACQ,"sw")
   SFO1 <- .varian.get_param(ACQ,"sfrq")
   SW <- SWH/SFO1
   ENDIAN = "big"

   # File Header: 32 bytes = 8x4
   to.read = file(FIDFILE,"rb")
   header<-readBin(to.read, what="integer",size=4, n=8, endian = ENDIAN)
   # Block Header: 28 bytes = 7x4
   block<-readBin(to.read, what="integer",size=4, n=7, endian = ENDIAN)
   ## Read FID
   signal<-readBin(to.read, what="integer",size=header[4], n=header[3], endian = ENDIAN)
   close(to.read)
   setwd(cur_dir)

   # Get the raw Real and Imaginary spectra 
   td <- length(signal)
   rawR <- signal[seq(from = 1, to = td, by = 2)]
   rawI <- signal[seq(from = 2, to = td, by = 2)]
   TD <- length(rawR)

   # null average
   mediar<-mean(as.integer(rawR[c((3*length(rawR)/4):length(rawR))]),na.rm = TRUE)
   mediai<--mean(as.integer(rawI[c((3*length(rawR)/4):length(rawR))]),na.rm = TRUE)
   rawR<-rawR-mediar
   rawI<-rawI-mediai
   fid <- rawR+1i*rawI

   GRPDLY <- 0
   PROBE <- .varian.get_param(ACQ,"probe_",type="string")
   SOLVENT <-  .varian.get_param(ACQ,"solvent",type="string")
   PULSE <- .varian.get_param(ACQ,"pslabel",type="string")
   ORIGPATH <- .varian.get_param(ACQ,"exppath",type="string")
   RELAXDELAY  <- .varian.get_param(ACQ,"d1")
   PULSEWIDTH <- .varian.get_param(ACQ,"pw90")
   SPINNINGRATE <- .varian.get_param(ACQ,"spin")
   TEMP <- .varian.get_param(ACQ,"temp") + 274.15

   proc <- list( INSTRUMENT="VARIAN", SOFTWARE="VnmrJ", ORIGIN="VARIAN", ORIGPATH=ORIGPATH, PROBE=PROBE, PULSE=PULSE, SOLVENT=SOLVENT, 
                RELAXDELAY=RELAXDELAY, SPINNINGRATE=SPINNINGRATE, PULSEWIDTH=PULSEWIDTH, TEMP=TEMP, 
                TD=TD, SW=SW, SWH=SWH, SFO1=SFO1, GRPDLY=GRPDLY, phc0=0, phc1=0 )
   spec <- list( path=DIR, proc=proc, fid=fid )

   spec

}

#### Read FID data and the main parameters needed to generate the real spectrum
#--  external routine
#-- filename: the nmrML file
.read.FID.nmrML <- function(filename)
{
   if (!file.exists(filename))
       stop("nmrML File ", filename, " does not exist\n")

   what <- "double"
   endian <- "little"
   sizeof <- 8
   compression <- "gzip"
   
   tree <- xmlTreeParse(filename)
   root <- xmlRoot(tree)
   fidData <- xmlElementsByTagName(root, "fidData", recursive = TRUE)[["acquisition.acquisition1D.fidData"]]
   b64string <- gsub("\n", "", xmlValue(fidData))
   byteFormat <- xmlAttrs(fidData)["byteFormat"]
   raws <- memDecompress(base64decode(b64string), type=compression)
   signal <- readBin(raws, n=length(raws), what=what, size=sizeof, endian = endian)
   TD <- length(signal)

   SFO1 <- as.double(xmlAttrs(xmlElementsByTagName(root, "irradiationFrequency", recursive = TRUE)[[1]])["value"])
   SWH <-  as.double(xmlAttrs(xmlElementsByTagName(root, "sweepWidth", recursive = TRUE)[[1]])["value"])
   SW <- SWH/SFO1
   TD  <-  as.integer(xmlAttrs(xmlElementsByTagName(root, "DirectDimensionParameterSet", recursive = TRUE)[[1]])["numberOfDataPoints"])

   TEMP <- as.double(xmlAttrs(xmlElementsByTagName(root, "sampleAcquisitionTemperature", recursive = TRUE)[[1]])["value"])
   RELAXDELAY <- as.double(xmlAttrs(xmlElementsByTagName(root, "relaxationDelay", recursive = TRUE)[[1]])["value"])
   SPINNINGRATE <- as.double(xmlAttrs(xmlElementsByTagName(root, "spinningRate", recursive = TRUE)[[1]])["value"])
   PULSEWIDTH <- as.double(xmlAttrs(xmlElementsByTagName(root, "pulseWidth", recursive = TRUE)[[1]])["value"])

   instrument <- xmlElementsByTagName(root, "instrumentConfiguration", recursive = TRUE)[[1]]
   INSTRUMENT <- xmlAttrs(xmlElementsByTagName(instrument,"cvParam")[[1]])["name"]
   PROBE <- xmlAttrs(xmlElementsByTagName(instrument,"userParam")[[1]])["value"]

   SOFTWARE <- ''
   PULSE   <- ''
   ORIGIN   <- '-'
   ORIGPATH <- '-'
   SOLVENT <- '-'

   td <- length(signal)
   rawR <- signal[seq(from = 1, to = td, by = 2)]
   rawI <- signal[seq(from = 2, to = td, by = 2)]

   mediar<-mean(as.integer(rawR[c((3*length(rawR)/4):length(rawR))]),na.rm = TRUE)
   mediai<--mean(as.integer(rawI[c((3*length(rawR)/4):length(rawR))]),na.rm = TRUE)
   rawR<-rawR-mediar
   rawI<-rawI-mediai
   fid <- rawR+1i*rawI

   # Fix the Group Delay parameter if Bruker instrument
   if ( regexpr("(B|b)ruker",INSTRUMENT)[[1]]>0 ) {
      ### Estimation of point number of the dead time for the FID (Group Delay)
      nd0 <-which.max(abs(Re(fid)))
      nd0p <- which.max(abs(Re(fid)[1:(nd0-1)]))
      if ((abs(Re(fid)[nd0p])/abs(Re(fid)[nd0]))>0.5) nd0 <- nd0p;
      while(sign(Re(fid)[nd0-1])==sign(Re(fid)[nd0])) nd0 <- nd0-1;
      nd0 <- nd0-1
      GRPDLY <- nd0
   } else {
      GRPDLY  <- 0
   }

   proc <- list( INSTRUMENT=INSTRUMENT, SOFTWARE=SOFTWARE, ORIGIN=ORIGIN, ORIGPATH=ORIGPATH, 
                PROBE=PROBE, PULSE=PULSE, SOLVENT=SOLVENT, TEMP=TEMP, 
                RELAXDELAY=RELAXDELAY, SPINNINGRATE=SPINNINGRATE, PULSEWIDTH=PULSEWIDTH,
                TD=TD, SW=SW, SWH=SWH, SFO1=SFO1, GRPDLY=GRPDLY, phc0=0, phc1=0 )
   spec <- list( path=filename, proc=proc, fid=fid )

   spec
}

#--  external routine
#-- DIR: bruker directory containing the FID file, 
.read.1r.bruker <- function(DIR, param=Spec1r.Procpar) {

   cur_dir <- getwd()
   setwd(DIR)
   
   # Acquisition parameters filename
   ACQFILE <- "acqus"
   if (!file.exists(ACQFILE)) 
       stop("Acquisition parameters File (acqus) does not exist\n")

   # Processing parameters filename
   SPECFILE <- paste(param$PDATA_DIR,"/1r",sep="")
   PROCFILE <- paste(param$PDATA_DIR,"/procs",sep="")
   if (!file.exists(PROCFILE)) 
       stop("Processing parameters File (procs) does not exist\n")

   # Read Acquisition parameters
   ACQ     <- readLines(ACQFILE)
   PROBE   <- .bruker.get_param(ACQ,"PROBHD",type="string")
   SOLVENT <- .bruker.get_param(ACQ,"SOLVENT",type="string")
   PULSE   <- .bruker.get_param(ACQ,"PULPROG",type="string")
   TEMP    <- .bruker.get_param(ACQ,"TE",type="string")
   RELAXDELAY  <- .bruker.get_param(ACQ,"D",type="string",  arrayType=TRUE)
   PULSEWIDTH <- .bruker.get_param(ACQ,"P",type="string",  arrayType=TRUE)
   SPINNINGRATE <- .bruker.get_param(ACQ,"MASR",type="string")
   TD      <- .bruker.get_param(ACQ,"TD")
   SW      <- .bruker.get_param(ACQ,"SW")
   SWH     <- .bruker.get_param(ACQ,"SW_h")
   SFO1    <- .bruker.get_param(ACQ,"SFO1")
   GRPDLY  <- .bruker.get_param(ACQ,"GRPDLY")
   DTYPA   <- .bruker.get_param(ACQ,"DTYPA")
   BYTORDA <- .bruker.get_param(ACQ,"BYTORDA")
   #ENDIAN = ifelse( BYTORDA==0, "little", "big")
   ENDIAN <- "little"
   INSTRUMENT <- "Bruker"
   SOFTWARE <- gsub("##TITLE= Parameter file, ","", gsub("\t\tVersion"," v.",ACQ[1]))
   ORIGIN   <- gsub("^[^=]+= ","", ACQ[which(simplify2array(regexpr("^..ORIGIN=",ACQ))>0)])
   ORIGPATH <- gsub("^.. ", "", ACQ[which(simplify2array(regexpr("acqus$",ACQ))>0)])
   SIZE = ifelse( DTYPA==0, 4, 8)
   GRPDLY <- tryCatch ( { if (!is.na(as.numeric(GRPDLY)) && as.numeric(GRPDLY)>0) round(GRPDLY); }, warning = function(w){}, error=function(e){} )
   if (is.null(GRPDLY)) {
        GRPDLY <- 0
   }

   # Read Processing parameters
   PROC <- readLines(PROCFILE)
   OFFSET <- .bruker.get_param(PROC,"OFFSET")
   SI <- .bruker.get_param(PROC,"SI")
   PHC0 <-  .bruker.get_param(PROC,"PHC0")
   PHC1 <-  .bruker.get_param(PROC,"PHC1")

   to.read = file(SPECFILE,"rb")
   signal<-rev(readBin(to.read, what="int",size=SIZE, n=SI, signed = TRUE, endian = ENDIAN))
   close(to.read)
   TD <- length(signal)

   setwd(cur_dir)

   dppm <- SW/(TD-1)
   pmax <- OFFSET
   pmin <- OFFSET - SW
   ppm <- seq(from=pmin, to=pmax, by=dppm)

   proc <- list( INSTRUMENT=INSTRUMENT, SOFTWARE=SOFTWARE, ORIGIN=ORIGIN, ORIGPATH=ORIGPATH, 
                PROBE=PROBE, PULSE=PULSE, SOLVENT=SOLVENT, TEMP=TEMP, 
                RELAXDELAY=RELAXDELAY, SPINNINGRATE=SPINNINGRATE, PULSEWIDTH=PULSEWIDTH,
                TD=TD, SW=SW, SWH=SWH, SFO1=SFO1, GRPDLY=GRPDLY, phc0=PHC0, phc1=PHC1, SI=SI )

   ### Fix the ppm range where the TSP signal could be present depending on the 'Spectrum Width'
   if (param$TSP) {
       param$TSP_PPM_RANGE <- 0.5*proc$SW - 4.5
   }

   proc$NegFID <- FALSE
   proc$NegPhi <- FALSE
   param$ZEROFILLING <- FALSE
   param$LINEBROADENING <- FALSE

   spec <- list( path=DIR, param=param, proc=proc, fid=NULL, pklist=NULL, I0=0, I1=0, int=signal, dppm=dppm, pmin=pmin, pmax=pmax, ppm=ppm )

   spec
}

#--------------------------------
# Pre-Processing
#--------------------------------

#### Apply some preprocessing: zero_filling, line broading
#--  Generate real and imaginary parts
.preprocess <- function(spec,param=Spec1r.Procpar)
{
    logfile <- param$LOGFILE
    ### Estimation of point number of the dead time for the FID (nd0)
    M0 <-which.max(abs(Re(spec$fid))); M<-M0;
    M0p <- which.max(abs(Re(spec$fid)[1:(M0-1)]))
    if ((abs(Re(spec$fid)[M0p])/abs(Re(spec$fid)[M0]))>0.5) M <- M0p;

    nd0 <- ifelse (spec$proc$GRPDLY>0, spec$proc$GRPDLY, 0)
    if (param$DEBUG) .v("GRPDLY = %d\n", nd0, logfile=logfile)
    
    ### Adjustment of the sign of the FID and of the phasing computation 
    ### so that the zero order phase (phc0) should be included in the angle range [0-180] degree
    spec$proc$NegFID <- FALSE
    spec$proc$NegPhi <- FALSE
    if (sign(Re(spec$fid)[M])<0) spec$proc$NegFID <- TRUE
    if (sign(Im(spec$fid)[which.max(abs(Im(spec$fid)))])==sign(Re(spec$fid)[M])) spec$proc$NegPhi <- TRUE
    if (spec$proc$NegFID) { signFID <- -1; } else { signFID <- 1; }

    ### shift FID of nd0 points
    if (nd0>0) {
       rawR<-signFID*c( Re(spec$fid)[-c(1:nd0)], rep(0,nd0) )
       rawI<-signFID*c( Im(spec$fid)[-c(1:nd0)], rep(0,nd0) )
    } else {
       rawR<-signFID*Re(spec$fid)
       rawI<-signFID*Im(spec$fid)
    }
    if (param$DEBUG) .v("TD = %d\n", length(rawR), logfile=logfile)

    ### TD needs to be power of 2; if not, apply a padding of zeros
    td <- length(rawR)
    tdp2 <- 2^round(log2(td)+0.4999)
    if (td < tdp2 ) {
       if(param$DEBUG) .v("Zero Padding = %d\n", tdp2 - td, logfile=logfile)
       rawR<-c( rawR, rep(0,(tdp2-td)) )
       rawI<-c( rawI, rep(0,(tdp2-td)) )
    }

    ### Zero filling
    if (param$ZEROFILLING) {
       TDMAX <- 131072
       #TDMAX <- 65536
       if(param$DEBUG) .v("Zero Filling (x%d)\n", round(TDMAX/length(rawR)), logfile=logfile)
       while ((length(rawR)/TDMAX) < 1 ) {
          td=length(rawR)
          rawR<-c( rawR, rep(0,td) )
          rawI<-c( rawI, rep(0,td) )
       }
    }
    td=length(rawR)
    if (param$DEBUG) .v("SI = %d\n", td, logfile=logfile)

    ### Line Broadening
    vlb <- 1
    if (param$LINEBROADENING) {
       if(param$DEBUG) .v("Line Broadening (LB=%f)\n", param$LB, logfile=logfile)
       vlb <- exp( -seq(0,td-1)*param$LB*pi/(2*spec$proc$SWH) )
    }
    fid<-vlb*(rawR+1i*rawI)
    fid[1] <- 0.5*fid[1]

    ### FFT of FID
    if(param$DEBUG) .v("FFT ...", logfile=logfile)
    fspec <- fft(fid)
    if(param$DEBUG) .v("OK\n", logfile=logfile)

    ### Rotation
    spec.p <- c( fspec[(td/2+1):td], fspec[1:(td/2)] )
    rawspec <- spec.p

    ### Fix the ppm range where the TSP signal could be present depending on the 'Spectrum Width'
    if (param$TSP) {
       param$TSP_PPM_RANGE <- 0.5*spec$proc$SW - 4.5
    }

    ### Save into the spec object instance
    spec$proc$TD <- length(rawspec)
    spec$param <- param
    spec$re <- Re(rawspec)
    spec$im <- Im(rawspec)
    spec$proc$phc0 <- 0
    spec$proc$phc1 <- 0
    spec$I0 <- 0
    spec$I1 <- 0
    spec$fid <- NULL

    ### return the spec object instance
    spec
}

.find_Peaks <- function(spec)
{
    logfile <- spec$param$LOGFILE
    if(spec$param$DEBUG) .v("Find Major Peaks ... ", logfile=logfile)

    spec.re <- 100*spec$int/max(abs(spec$int))

    # Init
    TD <- spec$proc$TD
    DW <- spec$proc$SW/(TD-1)
    FN=50
    minNratio0 <- spec$param$TSP_PPM_RANGE/spec$proc$SW
    minNratio1 <- ifelse( spec$proc$SW>15, 0.35, 0.3 )

    SIGMA <- 5
    N <- 101
    if (length(spec.re) < TD ) {
         spec.re <- c(spec.re, rep(0,TD-length(spec.re)))
    }
    TDC <- ifelse ( TD<=65536, TD-2, TD )
    V <- convolve(spec.re[1:TDC], SDL(seq(-N,N),SIGMA), conj=FALSE, type="open")[-c(1:(3*N))]
    sv <- simplify2array(lapply( c(3:61), function(x) { sd(V[((x-1)*TD/64):(x*TD/64)]); }))
    mv <- simplify2array(lapply( c(3:61), function(x) { mean(V[((x-1)*TD/64):(x*TD/64)]); }))
    sig<-min(sv)
    moy<- mv[which.min(sv)]
    V <- V - moy

    # Init. minInt0 & minInt (maxInt=100)
    minInt0 <- 0.5
    minInt <- 1
    # if TSP
    if (spec$param$TSP) {
        n <- round(spec$proc$TD*minNratio0)
        minInt0 <- round(max(abs(spec.re[1:n]))-0.001,4)
        minInt <- ifelse( minInt0<0.5, 2*minInt0, 1 )
    } else {
    # if PLASMA
       spec.tmp <- spec.re
       N <- which(spec.tmp==100)
       n <- round(abs(spec$proc$TD)*0.2/spec$proc$SW)
       spec.tmp[(N-n):(N+n)] <- 0
       minInt <- max(abs(spec.tmp))/100
       minInt0 <- 0.5*minInt
    }

    # Get Peak List
    pklist <- C_findPeaks(abs(spec.re), V, FN*sig, DW, minNratio0, minInt0, minInt, spec$param$PLASMA)

    # Find peaks for phasing
    I0 <- 1
    k <- 2;
    while ((pklist[[k]]$N/TD)<minNratio0) { 
          if ( pklist[[k]]$maxInt > pklist[[I0]]$maxInt ) { I0 <- k; }
          k <- k + 1;
    }
    while ((pklist[[k]]$N/TD)<minNratio1) { k <- k + 1; }
    I1 <- k
    for( i in (k+1):length(pklist) ) if (pklist[[i]]$maxInt>pklist[[I1]]$maxInt) { I1 <- i; }
    spec$I0 <- I0
    spec$I1 <- I1
    spec$pklist <- pklist
    if(spec$param$DEBUG) .v("OK\n", logfile=logfile)

    spec
}

.Optimphase0 <- function (spec)
{
   logfile <- spec$param$LOGFILE
   if (spec$param$DEBUG) .v("phc0 Optimization ... ", logfile=logfile)

   V<-90; X <- seq(-V,V,0.1);
   Y <- simplify2array(lapply( X, function(x) { C_phc0_fn(spec,x); }));
   Ys <- which(Y==0);
   spec$proc$phc0 <- ifelse( length(Ys)>0, 0.5*(X[Ys[length(Ys)]]+X[Ys[1]]), X[which.min(Y)] )
   V<-spec$proc$phc0; X <- seq(V-0.1,V+0.1,0.01);
   Y <- simplify2array(lapply( X, function(x) { C_phc0_fn(spec,x); }));
   Ys <- which(Y==0);
   spec$proc$phc0 <- ifelse( length(Ys)>0, 0.5*(X[Ys[length(Ys)]]+X[Ys[1]]), X[which.min(Y)] )

   if(spec$param$DEBUG) .v("OK: phc0 = %f\n",spec$proc$phc0, logfile=logfile)
   spec
}

.Optimphase1 <- function (spec)
{
   evalIM <- function (spec, x) {
             spec$phc1 <- x
             TD <- spec$proc$TD
             V <- C_corr_spec_im(spec)
             sv <- simplify2array(lapply( c(3:29), function(x) { sd(V[((x-1)*TD/32):(x*TD/32)]); }))
             mv <- simplify2array(lapply( c(3:29), function(x) { mean(V[((x-1)*TD/32):(x*TD/32)]); }))
             moy<- mv[which.min(sv)]
             V <- V - moy
             sum(V)
   }

   logfile <- spec$param$LOGFILE
   if (spec$param$DEBUG) .v("phc1 Optimization ... ", logfile=logfile)
   if (spec$param$PLASMA) {
       spec$proc$phc1 <- C_optim_phc(spec, c(-0.001, 0.001), C_phc1_fn)
   } else {
       #spec$proc$phc1 <- C_optim_phc(spec, c(-0.003, 0.003), C_phc1_fn)
       V<-0.008;
       X <- seq(-V,V,0.00001);
       # Y <- order(simplify2array(lapply( X, function(x) { C_phc1_fn(spec,x); })))[1:10];
       # spec$proc$phc1 <- X[ Y[ which.min(abs(X[Y])) ] ]
       W <- simplify2array(lapply( X, function(x) { C_phc1_fn(spec,x); }))
       Xm <- X[order(W)[1:4]][order(abs(X[order(W)[1:4]]))[1:2]]
       spec$proc$phc1 <- Xm[ which.min(abs(simplify2array(lapply( Xm, function(x) { evalIM(spec,x); })))) ]
   }

   if(spec$param$DEBUG) .v("OK: phc1 = %f\n",spec$proc$phc1, logfile=logfile)
   spec
}


### FID Processing - Main routine
#--    DIR        : absolute path of the Bruker/Varian directory
#--    procparams : list of processing parameters, see the  default parameter list 'Spec1rFromFID.params'
###
.CALL <- function ( Input, proc=Spec1r.Procpar )
{

   logfile <- proc$LOGFILE
   if(proc$DEBUG)
       if(proc$INPUT_SIGNAL == "fid") {
         .v("Read the FID ...",logfile=logfile)
       } else {
          .v("Read the 1R ...",logfile=logfile)
       }

   ## Read FID or 1r and parameters
   repeat {
      if (proc$VENDOR == "nmrml") {
         spec <- .read.FID.nmrML(Input)
         break
      }
      if (proc$VENDOR == "bruker" && proc$INPUT_SIGNAL == "fid") {
         spec <- .read.FID.bruker(Input)
         break
      }
      if (proc$VENDOR == "bruker" && proc$INPUT_SIGNAL == "1r") {
         spec <- .read.1r.bruker(Input,proc)
         break
      }
      if (proc$VENDOR == "varian"){
         proc$INPUT_SIGNAL <- "fid"
         spec <- .read.FID.varian(Input)
         break
      }
      break
   }

   repeat {
      if (proc$READ_RAW_ONLY) break

      if(proc$DEBUG) .v("OK\n",logfile=logfile)

      if ( proc$INPUT_SIGNAL == "fid") {
          ## Pre-processing: group delay, zero filling, line broadening
          spec <- .preprocess(spec,proc)

          ## Find Peaks
          spec$int <- spec$re
          spec <- .find_Peaks(spec)
          ## Phasing
          spec <- .Optimphase0(spec)
          spec <- .Optimphase1(spec)

          # Re-adjust peaks location
          spec.re <- C_corr_spec_re(spec)
          spec$pklist[[spec$I0]]$N <-  spec$pklist[[spec$I0]]$n1-1 + which.max(spec.re[(spec$pklist[[spec$I0]]$n1):(spec$pklist[[spec$I0]]$n2)])
          spec$pklist[[spec$I1]]$N <-  spec$pklist[[spec$I1]]$n1-1 + which.max(spec.re[(spec$pklist[[spec$I1]]$n1):(spec$pklist[[spec$I1]]$n2)])
          spec$int <- C_corr_spec_re(spec)
          spec$re <- NULL          
          spec$im <- NULL
      }

      if (proc$INPUT_SIGNAL == "1r") {
          ## Find Peaks
          spec <- .find_Peaks(spec)
          spec$pCAL <- (spec$pklist[[spec$I1]]$N-1)*spec$dppm + spec$pmin
          spec$Imax <- spec$int[spec$pklist[[spec$I1]]$N]
      }

      # Define spec list as a Spectrum object
      class(spec) = "Spectrum"
      # Ouput spec object

      break
   }

   flush.console()
   spec

}


### Write a Matrix of Spectrum in a binary mode (PACK format)
#   specMat : the Matrix of Spectrum : 1 row <=> 1 spectrum, 1 column <=> a same value of ppm
#   ppm_min, ppm_max :  the ppm range of the spectra
#   filepack : the full path of binary file
writeSpecMatrix = function(specMat, ppm_min, ppm_max, filepack)
{
   C_write_pack(specMat, ppm_min, ppm_max, filepack)
}

### Read a Matrix of Spectrum in a binary mode (PACK format)
#   Input: filepack : the full path of binary file
#   Output: a list with :
#         int: the Matrix of Spectrum : 1 row <=> 1 spectrum, 1 column <=> a same value of ppm
#         nspec & size : respectively the number of spectra and their size points
#         ppm_min & ppm_max : respectively the minimum and the maximum of the PPM range
readSpecMatrix = function(filepack)
{
   C_read_pack(filepack)
}
