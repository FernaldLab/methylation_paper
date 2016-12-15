# functions related to bsseq and BSmooth

################################################################################################
#### load methratio.py output files from multiple subjects
## Args
##  filename: generic filename for methratio.py output files
##  subjects: vector of names of subject directories that hold methratio.py output files
##  colClasses: vector defining colClasses argument to read.table()
##   for example, c('character','numeric',rep('NULL',3),'numeric','numeric',rep('NULL',5))
##  colNames: vector of names for columns read in from methratio.py output
##   for example, c('chr','pos','Cov','M');
##  header: logical indicating whether methratio.py output has header line
##  checkTime: logical indicating whether to print Sys.time() before each subject
## Output
##  list containing one data frame per subject
## Comment
##  Assumes working directory contains one sub-directory for each subject
##   each subject dir contains a methratio.py output file, named "filename" (arg 1)

.loadMethRatioData = function (filename, subjects, colClasses, colNames, header, checkTime=F, ...) {
  d = list();
  for (s in subjects) {
    if (checkTime) { print(Sys.time()) }
    cat(s,'\n');
    toRead = paste(s, '/', filename, sep='');
    d[[s]] = read.table(toRead, header=header, sep='\t', colClasses=colClasses, ...);
    names(d[[s]]) = colNames;
    gc();
  }
  return(d);
}

################################################################################################
#### filter based on coverage, then set coverage equal to methylation where coverage < methylation
## Args
##  datDF: data frame with coverage and methylation values, usually one element of output from .loadMethRatioData()
##   must have a column named "Cov" (coverage) and one named "M" (methylation)
##  reqCov: number indicating required coverage level
## Output
##  data frame that is filtered/corrected version of datDF
## Comment
##  "Coverage" is understood to mean the effective_CT count output from methratio.py
##   which is sometimes lower than the methylation value

.correctEffectiveCTAndFilterCoverage = function (datDF, reqCov) {
  if (!all(c('Cov','M') %in% names(datDF))) {
    stop('Input must have columns named "Cov" (for coverage) and "M" (for methylation)');
  }
  checkCov = datDF$Cov >= reqCov;
  datDF_f = datDF[checkCov, ];
  rows = which(datDF_f$M > datDF_f$Cov);
  datDF_f$Cov[rows] = datDF_f$M[rows];
  return(datDF_f);
}

################################################################################################
#### build a BSseq object from a data frame
## Args
##  df: data frame, usually one element of output from .loadMethRatioData() 
##   must at least contain columns named 'chr', 'pos', 'Cov', 'M'
##   can optionally also have a 'strand' column
##  reqCov: number indicating required coverage level
##  subject: name of subject, for sampleNames and pData elements of BSseq object
##  group: group of subject, for pData element of BSseq object
## Output
##  BSseq object
## Comment

.makeBSseqFromDF = function (df, reqCov, subject, group) {
  df_f = .correctEffectiveCTAndFilterCoverage(df, reqCov);
  if ('strand' %in% names(df_f)) {
    st = df_f$strand;
  } else {
    st = rep('*', nrow(df_f));
  }
  forGr = data.frame(chr=df_f$chr, start=df_f$pos, end=df_f$pos, strand=st);
  bsd = BSseq(M=as.matrix(df_f$M, ncol=1),
              Cov=as.matrix(df_f$Cov, ncol=1),
              gr=data.frame2GRanges(forGr),
              sampleNames=subject,
              pData=data.frame(group=group, row.names=subject)
  );
  return(bsd);
}

################################################################################################
#### make a list of BSseq objects from a list of data frames
## Args
##  dfList: list of data frames to pass to .makeBSseqFromDF(), usually output from .loadMethRatioData() 
##  reqCov: number indicating required coverage level
##  groups: named vector where values are group designations and names are subjects
##   names should be the same as names for dfList, will be used to build pData elements of BSseq objects
##  checkTime: logical indicating whether to print Sys.time() before each subject
## Output
##  list containing one BSseq object for every data frame in dfList
## Comment

.makeBSseqListFromDFList = function (dfList, reqCov, groups, checkTime=F) {
  bsdList = list();
  for (s in 1:length(dfList)) {
    if (checkTime) { print(Sys.time()) }
    this_dat = dfList[[s]];
    this_sub = names(dfList)[s]; 
    cat(this_sub,'\n');
    this_gp = groups[match(this_sub, names(groups))];
    bsdList[[this_sub]] = .makeBSseqFromDF(df=this_dat, reqCov=reqCov, subject=this_sub, group=this_gp);
  }
  return(bsdList);
}

################################################################################################
#### combine BSseq objects that are in a list into a single object
## Args
##  bsdList: list of BSseq objects, usually output from .makeBSseqListFromDFList()
##  checkTime: logical indicating whether to print Sys.time() before each subject
## Output
##  BSseq object 
## Comment

.combineBSDList = function (bsdList, checkTime=F) {
  if (checkTime) { print(Sys.time()) }
  nSubs = length(bsdList);
  bsdall = combine(bsdList[[1]], bsdList[[2]]);
  if (nSubs > 2) {
    if (checkTime) { print(Sys.time()) }
    for (s in 3:nSubs) {
      bsdall = combine(bsdall, bsdList[[s]]); gc();
    }
  }
  return(bsdall);
}

################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
##  

# .summaryBSseq = function (BSseqObj) {
# 
# }

################################################################################################
#### remove seqlevels that have been filtered out
## Args
##  BSD: BSseq object
## Output
##  BSseq object with 0-count seqlevels removed
## Comment
##

.fixSeqLevels = function (BSD) {
  seqnameCounts = table(seqnames(BSD));
  seqlevels(BSD) = names(seqnameCounts)[seqnameCounts>0];
  return(BSD);
}

################################################################################################
#### filter BSseq object to remove scaffolds before some number
####  seqnames must be of format "scaffold_1"
## Args
##  BSD: BSseq object
##  scaffoldNumber: number to filter beyond
## Output
##  filtered BSseq object
## Comment
##  very specific for Burtoni, wrote because .BSmoothAcrossScaffolds() often quits partway through
##  e.g. if .BSmoothAcrossScaffolds() quit on scaffold_100, would run:
##   BSDpost100 = .filterBeyondScaffold(BSD, 99)

.filterBeyondScaffold = function (BSD, scaffoldNumber) {
  if (!is.numeric(scaffoldNumber)) { stop(paste('Arg scaffoldNumber must be numeric', sep=''))  }
  seqsplit = unlist(strsplit(as.vector(seqnames(BSD)), '_'));
  seqsplitnums = suppressWarnings(as.numeric(seqsplit));
  seqsplitnums = seqsplitnums[!is.na(seqsplitnums)];
  filteredBSD = BSD[seqsplitnums > scaffoldNumber];
  return(.fixSeqLevels(filteredBSD));
}

################################################################################################
#### run BSmooth on all seqlevels in a BSseq object separately
## Args
##  BSD: BSseq object, usually output from using combine() to gather multiple samples into one object
##  NS, H, MAXGAP: parameters to BSmooth
## Output
##  nothing saved into workspace, writes one file per seqlevel to directory named using arg values 
##  seqlevels with fewer loci than NS will be smoothed with reduced NS value and saved in subdirectory '_modified_ns'
##  seqlevels with no loci will be skipped and written to a text file named '_noLoci.txt'
## Comment
##  RStudio sometimes hangs at around 3000 iterations, not clear why, be careful if leave running overnight, etc
##  add option to skip scaffolds that contain too few loci

.BSmoothAcrossScaffolds = function (BSD, NS, H, MAXGAP, skipLowNS=F, dirPrefix=NULL) {
  if (is.character(dirPrefix)) {
    DIR = paste(dirPrefix, '_ns', NS, 'h', H, 'mg', MAXGAP, sep='');
  } else {
    # get name of input BSseq object
    BSDname = deparse(substitute(BSD));
    # create output directory
    DIR = paste(BSDname, '_ns', NS, 'h', H, 'mg', MAXGAP, sep='');
  }
  suppressWarnings(dir.create(DIR));
  # if needed, make subdirectory for scaffolds with fewer loci than NS parameter
  if (any(table(seqnames(BSD)) < NS) & !skipLowNS) { suppressWarnings(dir.create(paste(DIR, '/_modified_ns', sep=''))) }
  skipped = c();
  # loop through scaffolds
  for (scaffold in seqlevels(BSD)) {
    # make BSseq object for current scaffold and fix seqlevels/seqnames
    tmpBSD = subset(BSD, seqnames(BSD)==scaffold);
    if (length(tmpBSD) == 0) {
      warning(paste('Skipped ', scaffold, ' because no methylation loci', sep=''));
      skipped = c(skipped, scaffold);
      next;
    }
    if (skipLowNS & (length(tmpBSD) < NS)) {
      warning(paste('Skipped ', scaffold, ' because too few methylation loci', sep=''));
      skipped = c(skipped, scaffold);
      next;
    }
    seqlevels(tmpBSD) = scaffold;
    # check if number of methylation loci on scaffold is >= ns parameter to BSmooth
    if (length(tmpBSD) >= NS) {
      numCs = NS;
      okNS = TRUE;
    } else {
      # reduce ns if needed
      numCs = length(tmpBSD);
      okNS = FALSE;
    }
    cat(paste('============================================================== ', scaffold, '\n', sep=''));
    # run BSmooth
    tmpSmooth = BSmooth(tmpBSD, ns=numCs, h=H, maxGap=MAXGAP, mc.cores=4, parallelBy='sample');
    # assign results to new variable named after scaffold
    outName = scaffold;
    assign(outName, tmpSmooth);
    # build output filename
    if (okNS) {
      filename = paste(DIR, '/', scaffold, '.RData', sep='');
    } else {
      filename = paste(DIR, '/_modified_ns/', scaffold, '_ns', numCs, '.RData', sep='');
    }
    save(list=outName, file=filename);
    gc()
  }
  write.table(skipped, file=paste(DIR, '/_noLoci.txt', sep=''), quote=F, col.names=F, row.names=F);
}

################################################################################################
#### 
## Args
##  dir: directory containing smoothed BSseq objects, usually output from .BSmoothAcrossScaffolds()
##  printEvery: number telling function when to print progress 
## Output
##  list containing a separate smoothed BSseq object for each seqlevel/scaffold
## Comment
##

.loadSmoothedScaffolds = function (dir, printEvery=NULL) {
  files = list.files(dir);
  files = files[grep('^scaffold', files)];
  fits = vector(mode='list', length(length(files)));
  if (is.numeric(printEvery)) {
    for (f in 1:length(files)) { 
      if (f %% printEvery == 0) { cat(paste(signif(f/length(files)*100,2), '%... ', sep='')) }
      fits[[f]] = get(load(paste(dir,'/',files[f],sep=''))) 
    }
  } else {
    for (f in 1:length(files)) { 
      fits[[f]] = get(load(paste(dir,'/',files[f],sep=''))) 
    }
  }
  names(fits) = gsub('.RData','',files);
  return(fits);
}

# ################################################################################################
# #### filter loci with certain level of coverage in every sample
# ## Args
# ##  FIT: smoothed BSseq object
# ##  reqCov: number representing required coverage level
# ##  printEvery: number telling function when to print progress 
# ## Output
# ##  returns version of FIT with ranges removed where >=1 subject had <reqCov coverage
# ## Comment
# ##
# 
# .filterLowCoverage = function (FIT, reqCov=5, printEvery=NULL) {
#   Cov = getCoverage(FIT);
#   keep = c();
#   if (is.numeric(printEvery)) {
#     for (row in 1:nrow(Cov)) {
#       if (row %% printEvery == 0) { cat(paste(signif(row/nrow(Cov)*100,2), '% ', sep='')) }
#       if (all(Cov[row, ] >= reqCov)) { keep = c(keep, row) }
#     }
#   } else {
#     for (row in 1:nrow(Cov)) {
#       if (all(Cov[row, ] >= reqCov)) { keep = c(keep, row) }
#     }
#   }
#   return(FIT[keep, ]);
# }

################################################################################################
#### filter loci with certain level of coverage in every sample
## Args
##  FIT: BSseq object
##  reqCov: number representing required coverage level
##  fix.seqlevels: logical indicating whether to remove seqlevels that had all loci filtered out
## Output
##  returns version of FIT with ranges removed where >=1 subject had <reqCov coverage
## Comment
##

.filterLowCoverage2 = function (FIT, reqCov=5, fix.seqlevels=T) {
  Cov = getCoverage(FIT);
  keep = apply(Cov, 1, function(f) all(f >= reqCov));
  FITnew = FIT[keep, ];
  if (fix.seqlevels) { FITnew = .fixSeqLevels(BSD=FITnew) }
  return(FITnew);
}

################################################################################################
#### filter a list of BSseq objects based on minimum required coverage, i.e. run .filterLowCoverage() across a list
#### checks whether each component has a certain number of loci remaining
## Args
##  FITS: list of smoothed BSseq objects, is assumed smoothing parameters were the same for all
##  reqCov: number representing required coverage level
##  ns.BSmooth: number representing value of BSmooth() ns parameter when data was smoothed 
## Output
##  list with 3 elements:
##   1- filtered version of FITS
##   2- vector of seqlevels/scaffolds that had all loci removed by filtering
##   3- vector of seqlevels/scaffolds that originally had >ns.BSmooth loci but now have <ns.BSmooth
## Comment
##  add option for function to pull ns.BSmooth from the data
##  fix progress reporting, which is overly verbose

.filterLowCoverageList = function (FITS, reqCov=5, ns.BSmooth=70) {
  filteredList = FITS;
  totalLoci = sum(sapply(FITS,length));
  noLociLeft = c();
  numLociBelowNS = c();
  cat('Working...\n');
  filteredList[[1]] = .filterLowCoverage2(filteredList[[1]], reqCov=reqCov);
  progress = 0;
  for (i in 2:length(filteredList)) {
    pr = signif(sum(sapply(FITS[1:(i-1)],length)) / totalLoci * 100, 2);
    if ((pr != progress) & (pr %% 5 == 0)) { 
      progress = pr;
      cat(paste(progress, '%... ', sep=''));
    }
    filteredList[[i]] = .filterLowCoverage2(filteredList[[i]], reqCov=reqCov, fix.seqlevels=F);
    if (length(filteredList[[i]]) < ns.BSmooth) {
      numLociBelowNS = c(numLociBelowNS, names(filteredList)[i]);
    }
    if (length(filteredList[[i]]) == 0) {
      noLociLeft = c(noLociLeft, names(filteredList)[i]);
      warning(paste('Filtering at ', reqCov, 'x removed all loci from ', names(filteredList)[i], sep=''));
    }
  }
  return(list(filtered=filteredList, noLociLeft=noLociLeft, numLociBelowNS=numLociBelowNS));
}

################################################################################################
#### get coefficient of variation (CV) for each column of methylation values in a BSseq object
## Args
##  FIT: BSseq object, usually smoothed but doesn't need to be
##  type: string indicating whether to extract "smooth" or "raw" methylation data from FIT
## Output
##  vector containing CV of samples in FIT
## Comment
##

.getMethylationCV = function (FIT, type='smooth') {
  methLvls = getMeth(FIT, type=type);
  return(apply(methLvls, 2, function(f) sd(f,na.rm=T)/mean(f,na.rm=T)));
}

################################################################################################
#### 
## Args
##  FIT: BSseq object, usually smoothed but doesn't need to be
##  type: string indicating whether to extract "smooth" or "raw" methylation data from FIT
## Output
##  vector containing CV of samples in FIT
## Comment
##

.getMethylationSummary = function (FIT, type='smooth') {
  methLvls = getMeth(FIT, type=type);
  return(apply(methLvls, 2, function(f) summary(f,na.rm=T)));
}

################################################################################################
#### get coefficient of variation (CV) for each column of methylation values in multiple BSseq objects
## Args
##  FITS: list of BSseq objects, all need to have exactly the same sampleNames
##  type: string indicating whether to extract "smooth" or "raw" methylation data from elements of FITS
## Output
##  data frame containing CV of samples in each element of FITS
##   rownames are names of FITS, colnames are sampleNames from elements of FITS
## Comment
##

.getMethylationCVList = function (FITS, type='smooth') {
  subjects = sampleNames(FITS[[1]]);
  CVs = as.data.frame(matrix(ncol=length(subjects), nrow=length(FITS)));
  dimnames(CVs) = list(names(FITS), subjects);
  for (i in 1:length(FITS)) {
    if (any(sampleNames(FITS[[i]]) != subjects)) { stop('All list elements must have same sampleNames') }
    CVs[i, ] = .getMethylationCV(FITS[[i]], type=type);
  }
  return(CVs);
}

################################################################################################
#### 
## Args
##  FITS: list of BSseq objects, all need to have exactly the same sampleNames
##  type: string indicating whether to extract "smooth" or "raw" methylation data from elements of FITS
## Output
##  data frame containing CV of samples in each element of FITS
##   rownames are names of FITS, colnames are sampleNames from elements of FITS
## Comment
##

.getMethylationSummaryList = function (FITS, type='smooth') {
  subjects = sampleNames(FITS[[1]]);
  statList = vector(mode='list', length=6);
  names(statList) = c('min','q1','median','mean','q3','max');
  for (i in 1:length(statList)) {
    statList[[i]] = as.data.frame(matrix(ncol=length(subjects), nrow=length(FITS)));
    dimnames(statList[[i]]) = list(names(FITS), subjects);
  }
  for (j in 1:length(FITS)) {
    if (any(sampleNames(FITS[[j]]) != subjects)) { stop('All list elements must have same sampleNames') }
    scafStats = as.data.frame(t(.getMethylationSummary(FITS[[j]], type=type)));
    statList$min[j, ] = scafStats$'Min.';
    statList$q1[j, ] = scafStats$'1st Qu.';
    statList$median[j, ] = scafStats$'Median';
    statList$mean[j, ] = scafStats$'Mean';
    statList$q3[j, ] = scafStats$'3rd Qu.';
    statList$max[j, ] = scafStats$'Max.';
  }
  return(statList);
}

################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
##

.getLociDistsForScaffolds = function (FITS, scaffoldLengths) {
  res = as.data.frame(matrix(ncol=6,nrow=length(FITS), dimnames=list(names(FITS),c('min','q1','med','mean','q3','max'))));
  res2 = res;
  for (i in 1:length(FITS)) {
    res[i, ] = summary(diff(as.data.frame(ranges(FITS[[i]]))[,1]));
    res2[i, ] = res[i, ] / scaffoldLengths[match(names(FITS)[i], scaffoldLengths[,1]), 2];
  }
  return(list(val=res, pct=res2));
}





################################################################################################
#### compute t-statistics across a list of smoothed, filtered BSseq objects
#### run BSmooth.tstat for each element in list
## Args
##  FITS.FiltList: list of smoothed and filtered BSseq objects
##  group1, group2: vectors of subjects in each group, must correspond to sampleNames of BSseq objects
##  estimate.var: variance estimate parameter to BSmooth.tstat()
##  printEvery: number telling function when to print progress 
## Output
##  list with 2 elements:
##   1- list of BSseqTstat objects corresponding to the objects in FITS.FiltList
##   2- vector of seqlevels/scaffolds with no loci
## Comment
##

.BSmooth.tstatList = function (FITS.FiltList, group1, group2, estimate.var, printEvery=NULL) {
  if (is.null(printEvery)) { printEvery = ceiling(length(FITS.FiltList)/10) }
  TList = FITS.FiltList;
  noData = c();
  cat('Working...\n');
  for (i in 1:length(TList)) {
    cat(paste0(names(TList)[i], '... '))
    if (is.numeric(printEvery)) {
     if (i %% printEvery == 5) { cat(paste(signif(i/length(TList)*100,2), '%... ', sep='')) }
    }
    if (length(TList[[i]]) == 0) {
      warning(paste('Skipping ', names(TList)[i], ' because no data, output $tstats[[i]] will be NA', sep=''));
      noData = c(noData, names(TList)[i]);
      TList[[i]] = NA;
      next;
    } 
    if (any(is.na(getMeth(TList[[i]])))) {
      warning(paste('Skipping ', names(TList)[i], ' because data includes NAs, output $tstats[[i]] will be NA', sep=''));
      noData = c(noData, names(TList)[i]);
      TList[[i]] = NA;
      next;
    }
    TList[[i]] = BSmooth.tstat(TList[[i]], group1=group1, group2=group2, estimate.var=estimate.var, 
                               local.correct=T, mc.cores=4, verbose=F
                               );
  }
  gc();
  return(list(tstats=TList, noData=noData));
}

################################################################################################
#### return statistics and ranges from BSseqTstat object in one dataframe
## Args
##  BSseqTstat: BSseqTstat object
## Output
##  dataframe containing ranges and statistics
## Comment
##

.getStatsFromBSseqTstat = function (BSseqTstat) {
  if (class(BSseqTstat) != 'BSseqTstat') {
    warning(paste(deparse(substitute(BSseqTstat)), ' is not a BSseqTstat object\n returning NA', sep=''));
    return(NA);
  }
  loci = as.data.frame(ranges(attributes(BSseqTstat)$gr));
  stats = as.data.frame(attributes(BSseqTstat)$stats);
  return(cbind(loci, stats));
}

################################################################################################
#### compute tstat cutoff values for each BSseqTstat object in a list
## Args
##  BSseqTstatList: list of BSseqTstat objects
##  stat: name of stat ('tstat' or 'tstat.corrected') to use for calculation
##  centile: desired percentile, e.g. 0.95 means you want 95th percentile values
## Output
##  dataframe where rows are scaffolds/seqlevels and columns are left and right 
##   critical values for desired percentile
## Comment
##  function gets both left and right cutoffs and warns about NAs

.getTstatQuantilesForEachScaffold = function (BSseqTstatList, stat='tstat.corrected', centile=0.95) {
  q1 = (1-centile)/2;
  cuts = c(q1, 1-q1);
  qlist = vector(mode='list', length=length(BSseqTstatList));
  names(qlist) = names(BSseqTstatList);
  for (i in 1:length(BSseqTstatList)) {
    stats = suppressWarnings(.getStatsFromBSseqTstat(BSseqTstatList[[i]]));
    if (!is.data.frame(stats)) { 
      warning(paste('Skipping ', names(qlist)[i], ', not a BSseqTstat object\n element, row ', i, ' of output will be NA', sep=''));
      qlist[[i]] = c(NA, NA);
      next;
    }
    statCol = match(stat, names(stats));
    if (any(is.na(stats[, statCol]))) {
      warning(paste('Removing NAs from ', names(qlist)[i], sep=''));
    }
    qlist[[i]] = quantile(stats[, statCol], cuts, na.rm=T);
  }
  return(as.data.frame(do.call('rbind',qlist)));
}

################################################################################################
#### compute tstat cutoff values for an entire list of BSseqTstat objects
## Args
##  BSseqTstatList: list of BSseqTstat objects
##  stat: name of stat ('tstat' or 'tstat.corrected') to use for calculation
##  centile: desired percentile, e.g. 0.95 means you want 95th percentile values
## Output
##  vector with left and right cutoff values 
## Comment
##  will ignore NAs and will not warn about any

.getTstatQuantilesAcrossScaffolds = function (BSseqTstatList, stat='tstat.corrected', centile=0.95) {
  q1 = (1-centile)/2;
  cuts = c(q1, 1-q1);
  ts = c();
  for (i in 1:length(BSseqTstatList)) {
    #print(names(BSseqTstatList)[i])
    if (class(BSseqTstatList[[i]]) != 'BSseqTstat') {
      warning(paste('Skipping ', names(BSseqTstatList)[i], ', is not a BSseqTstat object', sep=''));
      next;
    }
    stats = suppressWarnings(.getStatsFromBSseqTstat(BSseqTstatList[[i]]));
    statCol = match(stat, names(stats));
    ts = c(ts, stats[, statCol]);
  }
  return(quantile(ts, cuts, na.rm=T));
}

################################################################################################
#### find DMRs across a list of BSseqTstat objects, usually output from .BSmooth.tstatList()$tstats
## Args
##  TList: list of BSseqTstat objects
##  cutoff:
##  maxGap: maxGap parameter to dmrFinder()
##  stat: stat parameter to dmrFinder()
##  printEvery: number telling function when to print progress 
## Output
##  list with 2 elements:
##   1- list of dataframes with DMRs for each BSseqTstat object in TList
##   2- dataframe with all DMRs collapsed from 1
## Comment
##
.dmrFinderList = function (TList, cutoff, maxGap=300, stat='tstat.corrected', printEvery=NULL) {
  if (!is.numeric(cutoff)) { stop('Arg "cutoff" must be numeric') }
  if (all((length(cutoff)==1) & (cutoff < 1))) {
    sep_cuts = TRUE;
    cat(paste0('Computing ', cutoff, ' quantile t-stat cutoffs for each scaffold... \n might reduce power on smaller scaffolds...'));
    tstatcut = .getTstatQuantilesForEachScaffold(TList, stat=stat, centile=cutoff);
    if (any(rownames(tstatcut) != names(TList))) { stop('Something went wrong computing t-stat cutoffs') }
  } else if ((length(cutoff)==2)) {
    cat(paste0('Using ', cutoff[1], ' and ', cutoff[2], ' as t-stat cutoffs for all scaffolds\n'));
    sep_cuts = FALSE;
    tstatcut = cutoff;
  } else {
    stop();
  }
  if (is.null(printEvery)) { printEvery = floor(length(TList)/10) }
  dmrList = TList;
  numTests = length(dmrList);
  for (i in 1:numTests) {
    if (i %% printEvery == 5) { cat(paste(signif(i/length(dmrList)*100,2), '%... ', sep='')) }
    if (class(dmrList[[i]]) != 'BSseqTstat') {
      warning(paste(names(dmrList)[i], ' is not a BSseqTstat object, element ', i, ' of $dmrList will be NA', sep=''));
      dmrList[[i]] = NA;
      next;
    }
    stats = suppressWarnings(.getStatsFromBSseqTstat(dmrList[[i]]));
    thisStat = stats[, match(stat, names(stats))];
    ckNum = is.numeric(thisStat); ckNA = is.na(thisStat);
    if (any(!ckNum | ckNA)) {
      naRows = which(!ckNum | ckNA);
      if (length(naRows) == length(dmrList[[i]])) {
        warning(paste('Skipped ', names(dmrList)[i], ' because all rows had NAs', sep=''));
        dmrList[[i]] = NA;
        next;
      }
      dmrList[[i]] = dmrList[[i]][-naRows, ];
      warning(paste('Removed ', length(naRows), ' rows in ', names(dmrList)[i], ' due to NAs', sep=''));
    }
    if (sep_cuts) {
      these_cuts = as.numeric(tstatcut[i, ]);
      if (any(is.na(these_cuts))) {
        warning(paste('t-stat cutoffs are NA for ', names(dmrList)[i], sep=''));
        dmrList[[i]] = NA;
        next;
      } else {
        tmp = dmrFinder(dmrList[[i]], cutoff=these_cuts, maxGap=maxGap, stat=stat, verbose=F);
      }
    } else {
      tmp = dmrFinder(dmrList[[i]], cutoff=cutoff, maxGap=maxGap, stat=stat, verbose=F);
    }
    if (is.null(tmp)) {
      warning(paste('No DMRs for ', names(dmrList)[i], sep=''));
      dmrList[[i]] = NA;
    } else {
      dmrList[[i]] = tmp;
    }
  }
  dmrDF = do.call('rbind', dmrList); dmrDF = subset(dmrDF, !apply(dmrDF, 1, function(f) sum(is.na(f))==length(f)));
  return(list(dmrList=dmrList, dmrDF=as.data.frame(dmrDF), tstatcut=tstatcut));
}

################################################################################################
#### filter and order DMRs by column name/number
## Args
##  DMRs: dataframe of DMRs from dmrFinder()
##  filterOn: colname(s) or number(s) indicating columns to filter on
##  filterValues: numeric value(s) for filtering, must be same length as filterOn
##  greater: logical indicating whether filtering should remove values greater/less than filterValues
##  orderOn: colname or number indicating column to order by
##  decreasing: logical indicating whether ordering should arrange values in decreasing/increasing order
## Output
##  filtered, possibly ordered, DMR dataframe
## Comment
##  filtering only works on numeric values and always uses absolute values
##  ordering only works with one column and uses absolute values
##  if filtering removes all rows will print a warning and return NA
  
.filterAndOrderDMRs = function (DMRs, filterOn=NULL, filterValues=NULL, greater=T, orderOn=NULL, decreasing=T, verbose=T) {
  if (!is.character(filterOn) | !is.numeric(filterOn)) { filterOn = c('n', 'meanDiff') }
  if (!is.numeric(filterValues)) { filterValues = c(3, 0.1) }
  if (length(filterOn) != length(filterValues)) { stop('filterOn and filterValues must be same length') }
  if (is.character(filterOn)) {
    filterCols = match(filterOn, names(DMRs));     
  } else if (is.numeric(filterOn)) {
    filterCols = filterOn;
  } else {
    stop('filterOn must be character or numeric');
  }
  if (verbose) {
    testtype = ifelse(greater, ' at least ', ' no greater than ');
    cat(paste0('Filtering rows with', testtype, 
              paste(paste0(filterCols, '=', filterValues), 
                    collapse=', '), 
              ' (abs.vals)...\n'));
  }
  filteredDMRs = DMRs;
  for (fCol in 1:length(filterCols)) {
    if (greater) {
      keepRows = abs(filteredDMRs[,filterCols[fCol]]) >= filterValues[fCol];
    } else {
      keepRows = abs(filteredDMRs[,filterCols[fCol]]) <= filterValues[fCol];
    }
    if (sum(keepRows) == 0) { 
      warning(paste0('filtering on ', 
                     names(filteredDMRs)[filterCols[fCol]], '=', filterValues[fCol], 
                     ' removes all rows, returning NA'));
      return(NA);
    } else {
      filteredDMRs = filteredDMRs[keepRows, ];
    }
  }
  
  if (!is.null(orderOn)) {
    if (is.character(orderOn)) {
      orderCol = match(orderOn, names(filteredDMRs));
    } else if (is.numeric(orderOn)) {
      orderCol = orderOn;
    } else {
      stop('orderOn must be character or numeric');
    }
    filteredDMRs = filteredDMRs[order(abs(filteredDMRs[,orderCol]), decreasing=decreasing), ];
  }
  return(filteredDMRs);
}
  
################################################################################################
#### filter and order a list of DMRs by column name/number
## Args
##  DMRsList: list of DMR dataframes output by dmrFinderList()
##  filterOn: colname(s) or number(s) indicating columns to filter on
##  filterValues: numeric value(s) for filtering, must be same length as filterOn
##  greater: logical indicating whether filtering should remove values greater/less than filterValues
##  orderOn: colname or number indicating column to order by
##  decreasing: logical indicating whether ordering should arrange values in decreasing/increasing order
## Output
##  list with 2 elements:
##   1- list of filtered, possibly ordered, DMR dataframes
##   2- dataframe with all DMRs collapsed from 1, reordered if required by orderOn
## Comment
##  if filtering removes all rows in a given DMR dataframe that element of $filteredList will be NA
##   warnings will be printed when function finishes

.filterAndOrderDMRsList = function (DMRsList, filterOn=NULL, filterValues=NULL, greater=T, orderOn=NULL, decreasing=T, verbose=T) {
  if (!is.character(filterOn) | !is.numeric(filterOn)) { filterOn = c('n', 'meanDiff') }
  if (!is.numeric(filterValues)) { filterValues = c(3, 0.1) }
  if (length(filterOn) != length(filterValues)) { stop('filterOn and filterValues must be same length') }
  if (verbose) {
    if (is.character(filterOn)) {
      vCols = filterOn;
    } else if (is.numeric(filterOn)) {
      vCols = match(filterOn, names(DMRsList[[1]]));
    } else {
      stop('filterOn must be character or numeric');
    }
    testtype = ifelse(greater, ' at least ', ' no greater than ');
    cat(paste('Filtering rows with', testtype, 
              paste(paste(vCols, '=', filterValues, sep=''), collapse=', '), 
              ' (abs.vals)...\n', sep=''
              )
        );
  }
  newList = DMRsList;
  for (i in 1:length(newList)) {
    if (!is.data.frame(newList[[i]])) {
      warning(paste('Skipping ', names(newList)[i], ' because no data', sep=''));
      next;
    }
   # print(head(newList[[i]]))
    newList[[i]] = .filterAndOrderDMRs(newList[[i]], filterOn=filterOn, filterValues=filterValues, 
                                       greater=greater, orderOn=orderOn, decreasing=decreasing, verbose=F
                                       );
  }
  dmrDF = do.call('rbind', newList);
  dmrDF = subset(dmrDF, !apply(dmrDF, 1, function(f) sum(is.na(f))==length(f)));
  if (!is.null(orderOn)) {
    if (is.character(orderOn)) {
      orderCol = match(orderOn, names(dmrDF));
    } else if (is.numeric(orderOn)) {
      orderCol = orderOn;
    } else {
      stop('orderOn must be character or numeric');
    }
    dmrDF = dmrDF[order(abs(dmrDF[,orderCol]), decreasing=decreasing), ];
  }
  return(list(filteredList=newList, filteredListCollapsed=as.data.frame(dmrDF)));
}

################################################################################################
#### sort DMRs based on average rank across multiple columns
## Args
##  DMRs: dataframe of DMRs from dmrFinder()
##  rankCols: vector of column names to use in ranking
## Output
##  list with 2 elements:
##   1- dataframe of ranks, including column with mean rank
##   2- re-order DMRs
## Comment
##

.orderDMRsByCompositeRank = function (DMRs, rankCols=c('n','width','areaStat','maxStat','meanDiff')) {
  r = as.data.frame(matrix(nrow=nrow(DMRs), ncol=length(rankCols), dimnames=list(rownames(DMRs),rankCols)));
  dmrCols = match(rankCols, names(DMRs));
  rCol = 1;
  for (i in dmrCols) {
    r[,rCol] = rank(-abs(DMRs[,i]));
    rCol = rCol+1;
  }
  ranks = cbind(r, avgRank=apply(r, 1, mean));
  newOrder = order(ranks$avgRank);
  return(list(ranks=ranks[newOrder, ], rankedDMRs=DMRs[newOrder, ]));
}

################################################################################################
#### call plotManyRegions() to plot DMRs
## Args
##  DMRs: dataframe of DMRs from dmrFinder()
##  BSD: smoothed BSseq object that was input to BSmooth.tstat
##  BSD.Tstat: output of BSmooth.tstat, not required but tstats will be added to plots if provided
##  colors: vector of colors to represent smoothed values for different subjects, order corresponds to order of sampleNames(BSD)
##  extend: number representing how far out to plot on either side of DMR
##  addPoints: logical indicating whether to plot points representing actual methylation values
##  filename: name for output pdf file, if NULL filename will be set with names of DMRs and BSD
## Output
##  nothing saved to workspace unless DMRs is not a dataframe (see comment)
##  writes a pdf file containing the plots
## Comment
##  if DMRs is not a dataframe a warning will be printed and NA returned

.plotDMRs = function(DMRs, BSD, BSD.Tstat=NULL, colors=c('blue','red','red','blue'), extend=5000, addPoints=T, filename=NULL, wh=c(12,6), ...) {
  if (!is.data.frame(DMRs)) {
    warning('Returning NA because DMRs is not a dataframe');
    return(NA)
  }
  pData = pData(BSD);
  pData$col = colors;
  pData(BSD) = pData;
  print(pData(BSD));
  if (!is.null(filename)) {
    if (!is.character(filename)) { stop('filename must be character string') }
  } else {
    filename = paste(deparse(substitute(BSD)), '_', deparse(substitute(DMRs)), '.pdf', sep='');
  }
  pdf(file=filename, width=wh[1], height=wh[2]);
  plotManyRegions(BSD, DMRs, extend=extend, addRegions=DMRs, addPoints=addPoints, BSseqTstat=BSD.Tstat, ...);
  dev.off();
}


.plotSingleDMR = function (dmrs, row, fitlist, bpleft=0, bpright=0, colors=c('midnightblue','greenyellow','greenyellow','midnightblue'), ...) {
  dmr = dmrs[row, 1:16];print(bpright)
  chr = dmr$chr; dmr$start = dmr$start-bpleft; dmr$end = dmr$end+bpright;
  thisfit = fitlist[[which(names(fitlist) == chr)]];
  pData(thisfit)$col = colors;
  plotRegion(BSseq=thisfit, region=dmr, ...)
}



################################################################################################
#### call plotManyRegions() on each element of a list containing DMR dataframes
## Args
##  DMRsList: list of DMR dataframes, probably from .dmrFinderList() or .filterAndOrderDMRsList()
##  BSDList: list of smoothed BSseq objects that were input to .BSmooth.tstatList()
##  BSD.TstatList: list of BSseqTstat objects from .BSmooth.tstatList()
##  colors: vector of colors to represent smoothed values for different subjects, order corresponds to order of sampleNames in each list element
##  extend: number representing how far out to plot on either side of DMR
##  addPoints: logical indicating whether to plot points representing actual methylation values
##  DIR: name of directory to store pdf files, if NULL will be set to name of DMRsList
## Output
##  nothing saved to workspace, creates directory to hold one pdf file per list element
##   pdfs will be named with names of DMRsList elements
## Comment
##  names of DMRsList, BSDList, and BSD.TstatList (if provided) must be exactly the same or will throw error
##  will print warning for any element of DMRsList that is not a dataframe

.plotDMRsAcrossList = function (DMRsList, BSDList, BSD.TstatList=NULL, colors=c('blue','red','red','blue'), extend=5000, addPoints=T, DIR=NULL, ...) {
  if (any(names(DMRsList) != names(BSDList))) { 
    stop(paste('names of DMRsList and BSDList must be exactly the same', sep=''))
  }
  if (!is.null(BSD.TstatList)) {
    if (any(names(DMRsList) != names(BSD.TstatList))) { 
      stop(paste('names of DMRsList and BSD.TstatList must be exactly the same', sep='')) 
    }
  }
  if (is.null(DIR)) { DIR = deparse(substitute(DMRsList)) }
  dir.create(DIR);
  if (!is.null(BSD.TstatList)) {
    for (i in 1:length(DMRsList)) {
      if (!is.data.frame(DMRsList[[i]])) { next }
      cat(names(DMRsList)[i],'\n')
      .plotDMRs(DMRs=DMRsList[[i]], BSD=BSDList[[i]], BSD.Tstat=BSD.TstatList[[i]], 
                colors=colors, extend=extend, addPoints=addPoints, filename=paste(DIR, '/', names(DMRsList)[i], '.pdf', sep=''), ...
                )
    }
  } else {
    for (i in 1:length(DMRsList)) {
      if (!is.data.frame(DMRsList[[i]])) { next }
      cat(names(DMRsList)[i],'\n')
      .plotDMRs(DMRs=DMRsList[[i]], BSD=BSDList[[i]], BSD.Tstat=NULL, 
                colors=colors, extend=extend, addPoints=addPoints, filename=paste(DIR, '/', names(DMRsList)[i], '.pdf', sep=''), ...
                )
    }
  }
}

################################################################################################
#### call getMeth() on a specific range within a BSseq object 
## Args
##  chr: string representing seqlevel/scaffold that range is on
##  start: numeric starting position of range
##  end: numeric ending position of range
##  BSseq: BSseq object, probably smoothed
##  ... additional parameters to getMeth(), e.g. type='raw'
## Output
##  matrix of methylation values, rows are loci, columns are subjects
## Comment
##  by default tries to get smoothed methylation values

.getMethForDMRRange = function (chr, start, end, BSseq, ...) {
  rangeGR = GRanges(seqnames=chr, ranges=IRanges(start=start, end=end));
  return(getMeth(BSseq=BSseq, regions=rangeGR, ...)[[1]]);
}

################################################################################################
#### look for DMRs that seem to be driven by only one animal within a group 
## Args
##  DMRs: dataframe of DMRs from dmrFinder()
##  smoothedList: BSseq object, usually smoothed, if not must add arg type='raw'
##  method: string indicating type of test to do, must start with either 'd' or 'c'
##   if 'd', will compare differences between average methylation values across a DMR 
##   if 'c', will compare correlation of methylation values across a DMR
##  strict: logical indicating strictness of test
##   if TRUE, will test whether all within group diffs are less than all across group diffs
##   if FALSE, will test whether within group diffs are less than the average across group diff
##    if method='c', will test cors instead of diffs, where within group values should be higher
##  ... additional options to .getMethForDMRRange(), i.e. getMeth()
## Output
##  integer dataframe with 2 columns
##   all- each row is either 0 or 1, reflecting whether all tests returned the same value (1) or not (0)
##   count- each row value is number of comparisons that returned TRUE
## Comment
##  assumes pData for all Bsseq objects in smoothedList is the same
##  method='cor', strict=F, removes DMRs that it probably shouldn't

.checkMethCoherence = function (DMRs, smoothedList, method=c('diff','cor'), strict=F, ...) {
  if (!all(unique(DMRs$chr) %in% names(smoothedList))) { stop() }
  if (!any(grepl('^d|c', method))) { stop() }
  if (!is.logical(strict)) { stop() }
  # get group info 
  grp = pData(smoothedList[[1]])$group;         #print(grp)
  dmrs = DMRs[, match(c('chr','start','end'), names(DMRs))];
  
  #res = vector(length=nrow(dmrs));
  res = as.data.frame(matrix(nrow=nrow(DMRs), ncol=2, 
                             dimnames=list(rownames(DMRs), 
                                           c('all','count'))));
  for (row in 1:nrow(dmrs)) {   #print('===========================')
    range = dmrs[row, ];
    ind = match(range$chr, names(smoothedList));
    meth = .getMethForDMRRange(chr=range$chr, start=range$start, end=range$end, 
                               BSseq=smoothedList[[ind]], ...);
    
    if (grepl('^d', method)) {
      methAvgs = apply(meth, 2, mean);
      diffMat = outer(methAvgs, methAvgs, '-');        #print(diffMat)
      theseRes = c();
      for (i in 1:nrow(diffMat)) { 
        grpInds = which(grp == grp[i]);
        othInds = which(grp != grp[i]);
        grpDiff = abs(diffMat[i, grpInds[grpInds!=i]]);  #print(grpDiff)
        othDiffs = abs(diffMat[i, othInds]);            # print(othDiffs)
        if (strict) {
          theseRes = c(theseRes, ifelse(all(grpDiff < othDiffs), TRUE, FALSE));
        } else {
          theseRes = c(theseRes, ifelse(grpDiff < mean(othDiffs), TRUE, FALSE));
        } #; print(theseRes)
      }
    } else {
      corMat = abs(cor(meth));      #print(corMat)
      theseRes = c();
      for (i in 1:nrow(corMat)) {         #print('-----------------')
        grpInds = which(grp == grp[i]);
        othInds = which(grp != grp[i]);
        grpCor = corMat[i, grpInds[grpInds!=i]];   #print(grpCor)
        othCors = corMat[i, othInds];              #print(mean(othCors))
        if (strict) {
          theseRes = c(theseRes, ifelse(all(grpCor > othCors), TRUE, FALSE));
        } else {
          theseRes = c(theseRes, ifelse(grpCor > mean(othCors), TRUE, FALSE));
        } 
      }
    }
   # res[row] = ifelse(count, sum(theseRes), all(theseRes));
   # res[row] = all(theseRes);
    res[row, ] = c(all(theseRes), sum(theseRes));
  }
  return(res);
}

################################################################################################
#### compare two dataframes of DMRs to find exact matches based on chr, start, end columns
## Args
##  dmrs1, dmrs2: dataframes of DMRs from dmrFinder()
## Output
##  list with 3 elements:
##   1- list containing 2 dataframes, one holding matching DMRs from each input dataframe
##   2- filtered version of dmrs1 containing DMRs that don't match any in dmrs2
##   3- filtered version of dmrs2 containing DMRs that don't match any in dmrs1
## Comment
##  names of dmrs1 and dmrs2 must be identical and have columns named 'chr', 'start', and 'end'

.compareDMRsExact = function (dmrs1, dmrs2) {
  if (any(names(dmrs1) != names(dmrs2))) { stop('Names of dmrs1 and dmrs2 must exactly match') }
  nameCols = match(c('chr','start','end'), names(dmrs1));
  chr = nameCols[1]; start = nameCols[2]; end = nameCols[3];
  names1 = paste(paste(dmrs1[,chr], dmrs1[,start],sep=':'), dmrs1[,end], sep='-');
  names2 = paste(paste(dmrs2[,chr], dmrs2[,start],sep=':'), dmrs2[,end], sep='-');
  inCommon = intersect(names1, names2);
  in1not2 = setdiff(names1, names2);
  in2not1 = setdiff(names2, names1);
  common1 = dmrs1[names1 %in% inCommon, ];
  common2 = dmrs2[names2 %in% inCommon, ];
  only1 = dmrs1[names1 %in% in1not2, ];
  only2 = dmrs2[names2 %in% in2not1, ];
  return(list(both=list(common1=common1, common2=common2), only1=only1, only2=only2));
}


.compareRangesExact = function (r1, r2, type=c('gr','df')) {
  #if (any(names(dmrs1) != names(dmrs2))) { stop('Names of dmrs1 and dmrs2 must exactly match') }
  nameCols = match(c('chr','start','end'), names(dmrs1));
  chr = nameCols[1]; start = nameCols[2]; end = nameCols[3];
  names1 = paste(paste(dmrs1[,chr], dmrs1[,start],sep=':'), dmrs1[,end], sep='-');
  names2 = paste(paste(dmrs2[,chr], dmrs2[,start],sep=':'), dmrs2[,end], sep='-');
  inCommon = intersect(names1, names2);
  in1not2 = setdiff(names1, names2);
  in2not1 = setdiff(names2, names1);
  common1 = dmrs1[names1 %in% inCommon, ];
  common2 = dmrs2[names2 %in% inCommon, ];
  only1 = dmrs1[names1 %in% in1not2, ];
  only2 = dmrs2[names2 %in% in2not1, ];
  return(list(both=list(common1=common1, common2=common2), only1=only1, only2=only2));
}

################################################################################################
#### compare two dataframes of DMRs 
## Args
##  dmrs1, dmrs2: dataframes of DMRs from dmrFinder()
## Output
##  list with 3 elements:
##   1- list of 2 row dataframes that contain overlapping DMRs from dmrs1 and dmrs2
##   2- filtered version of dmrs1 containing DMRs that don't overlap with anything in dmrs2
##   3- filtered version of dmrs2 containing DMRs that don't overlap with anything in dmrs1
## Comment
##  names of dmrs1 and dmrs2 must be identical and have columns named 'chr', 'start', and 'end'

.compareDMRsOverlap = function (dmrs1, dmrs2, ...) {
  if (any(names(dmrs1) != names(dmrs2))) { stop('Names of dmrs1 and dmrs2 must exactly match') }
  if (!all(c('chr','start','end') %in% names(dmrs1))) { stop('Names must include "chr", "start", and "end"') }
  gr1 = GRanges(seqnames=dmrs1$chr, 
                ranges=IRanges(start=dmrs1$start, end=dmrs1$end));#print(gr1)
  gr2 = GRanges(seqnames=dmrs2$chr, 
                ranges=IRanges(start=dmrs2$start, end=dmrs2$end));
  overlaps1 = as.data.frame(suppressWarnings(findOverlaps(gr1, gr2, ...)));#print(overlaps1)
  overlaps2 = as.data.frame(suppressWarnings(findOverlaps(gr2, gr1, ...)));  ### Probably don't need second findOverlaps   
  ck1 = all(overlaps1[,1] == sort(overlaps2[,2]));                      ###  or to do these checks
  ck2 = all(overlaps2[,1] == sort(overlaps1[,2]));                      ###  but will be safe for now
  if (!ck1 | !ck2) {
    stop('findOverlap calls were not symmetric, something is unusual about input datasets');
  } else {
    only1 = setdiff(1:nrow(dmrs1), overlaps1$queryHits);
    only2 = setdiff(1:nrow(dmrs2), overlaps2$queryHits);
    ovlist = vector(mode='list', length=3);
    names(ovlist) = c('overlaps', 'only1', 'only2');
    for (i in 1:nrow(overlaps1)) {
      ovlist$overlaps[[i]] = rbind(dmrs1[overlaps1[i,1], ], dmrs2[overlaps1[i,2], ]);
    }
    ovlist$only1 = dmrs1[only1, ];
    ovlist$only2 = dmrs2[only2, ];
  }
  return(ovlist);
}

###


###



# # slow slow slow
# # should just fix this issue in .compareDMRsOverlap when building ovlist$overlaps
.getMultipleOverlapInds = function (ovlist, printEvery=100) {
  uints = unique(unlist(lapply(ovlist, function(f) paste0(f$chr,':',f$start,'-',f$end))));
  dinds = list();
  for (i in 1:length(uints)) {
    if (i %% printEvery == 0) { cat(i,'...') } 
    ind = which(sapply(ovlist, function(f) any(paste0(f$chr,':',f$start,'-',f$end) == uints[i])));
    if (length(ind) > 1) {
      dinds[[length(dinds)+1]] = ind;
    }
  }
  return(dinds);
}


.cleanMultipleOverlaps = function (ovlist, dinds=NULL) {
  if (is.null(dinds)) {
    cat('Finding DMRs with multiple overlaps...');
    dinds = .getMultipleOverlapInds(ovlist, printEvery=100);
  }
  nomult = 1:length(ovlist);
  nomult = setdiff(nomult, unlist(dinds));
  newovlist = ovlist[nomult];
  for (ivec in dinds) {
    thism = do.call('rbind', ovlist[ivec]);
    drows = duplicated(paste0(thism$chr,':',thism$start,'-',thism$end));
    newovlist[[length(newovlist)+1]] = thism[!drows, ]; 
  }
  return(newovlist);
}

################################################################################################
#### combine 2 DMR dataframes after calling .compareDMRsOverlap()
## Args
##  dmrs1, dmrs2: dataframes of DMRs from dmrFinder()
## Output
##  dataframe with combination of dmrs1 and dmrs2
## Comment
##  in cases where ranges in dmrs1 and dmrs2 overlap the row from dmrs1 will be used in the output
##   rest of output rows are ranges unique to dmrs1 or dmrs2

.combineOverlappedDMRs = function (dmrs1, dmrs2, dmrsOverlap=NULL) {
  ov0 = .compareDMRsOverlap(dmrs1, dmrs2);
  ov = do.call('rbind', lapply(ov0$overlaps, function(f) f[1,]));
  ov = rbind(ov, ov0$only1, ov0$only2);
  return(ov);
}


# dmr stats end up being from whatever dmr table was first in the list
# .getCommonDmrsAcrossDmrListExact = function (dmrlist) {
#   commonDmrs = .compareDMRsExact(dmrlist[[1]], dmrlist[[2]])$both[[1]];
#   for (d in 3:length(dmrlist)) {
#     commonDmrs = .compareDMRsExact(commonDmrs, dmrlist[[d]])$both[[1]];
#   }
#   return(commonDmrs);
# }

# assumes rows of dmrdf are known overlapping dmr intervals
#  i.e. scaffold number and direction are the same
.mergeOverlappingDmrs = function (dmrdf, returnAll=FALSE) {
  out = rbind(dmrdf, rep(NA, ncol(dmrdf)));
  mm = range(as.integer(c(dmrdf$start, dmrdf$end)));
  out[nrow(out), 1] = dmrdf$chr[1];
  out[nrow(out), 2:3] = c(min(mm), max(mm));
  out[nrow(out), 7:15] = apply(dmrdf[,7:15], 2, mean);
  out[nrow(out), 16] = dmrdf$direction[1];
  ind = rep(returnAll, nrow(out)-1);
  return(out[c(ind, TRUE), ]);
}

.mergeOverlappingDmrsList = function (dmrOvList, returnAll=FALSE) {
  tmpdf = .mergeOverlappingDmrs(dmrOvList[[1]], returnAll=returnAll);
  for (i in 2:length(dmrOvList)) {
    tmpdf = rbind(tmpdf, .mergeOverlappingDmrs(dmrOvList[[i]], returnAll=returnAll));
  }
  return(tmpdf);
}

.getCommonDmrsAcrossDmrListExact = function (dmrlist) {
  cat(paste0('Getting first set of matches...'));
  tmp = .compareDMRsExact(dmrlist[[1]], dmrlist[[2]])$both;
  df0 = as.data.frame(rbind(tmp$common1, tmp$common2));
  for (d in 3:length(dmrlist)) {
    cat(paste0(d, '...'));
    dmrints = paste0(df0$chr, df0$start, df0$end);
    tmp = .compareDMRsExact(df0[!duplicated(dmrints), ], dmrlist[[d]])$both$common2;
    tmpints = paste0(tmp$chr, tmp$start, tmp$end);
    df0 = as.data.frame(rbind(df0[dmrints %in% tmpints, ], tmp));
  }
  dmrints = paste0(df0$chr, ':', df0$start, '-', df0$end);
  cat('Merging...');
  dfsplit = split(df0, dmrints);
  if (any(sapply(dfsplit, nrow) != length(dmrlist))) { stop('Some DMRs not in all lists') }
  dfmerge = .mergeOverlappingDmrsList(dfsplit, returnAll=FALSE);
  outlist = list(all=df0[order(df0$chr, df0$start, df0$end), ],
                 merge=dfmerge[order(-abs(dfmerge$areaStat)), ]);
  return(outlist);
}

.getCommonDmrsAcrossDmrListOverlap = function (dmrlist, minoverlap=10L) {
  cat(paste0('Getting first set of overlaps...'));
  ovDmrs = .compareDMRsOverlap(dmrlist[[1]], dmrlist[[2]], minoverlap=minoverlap)$overlaps;
  ovDmrs = .cleanMultipleOverlaps(ovDmrs, dinds=NULL);
  ovDmrs = .mergeOverlappingDmrsList(ovDmrs, returnAll=FALSE);
  for (d in 3:length(dmrlist)) {
    cat(paste0(d, '...'));
    ovDmrs = .compareDMRsOverlap(ovDmrs, dmrlist[[d]], minoverlap=minoverlap)$overlaps;
    ovDmrs = .cleanMultipleOverlaps(ovDmrs, dinds=NULL);
    ovDmrs = .mergeOverlappingDmrsList(ovDmrs, returnAll=FALSE);
  }
  return(ovDmrs[order(-abs(ovDmrs$areaStat)), ]);
}


################################################################################################
#### draw box-and-whisker plots comparing specified columns from 2 sets of DMRs
## Args
##  dmrs1, dmrs2: dataframes of DMRs from dmrFinder()
##  colNames: vector of column names to compare
##  mfrow: 2 number vector to define c(row, col) of plotting grid
##  abs: logical indicating whether to convert all numbers to absolute values
##  ...  additional options to boxplot
## Output
##  nothing saved to workspace, draws to standard Quartz window
## Comment
##

.compareDMRsPlots = function (dmrs1, dmrs2, colNames=NULL, mfrow=NULL, abs=F, grpNames=NULL, ...) {
  if (!all(names(dmrs1) == names(dmrs2))) { stop('Names of dmrs1 and dmrs2 must match exactly') }
  if (!is.character(colNames)) {
    colNames = c('n','width','invdensity','areaStat','maxStat','meanDiff','group1.mean','group2.mean','tstat.sd');
  }
  if (!is.numeric(mfrow) | length(mfrow)!=2) { stop('mfrow must be a 2 number vector') }
  if (is.null(grpNames)) { grpNames = c(deparse(substitute(dmrs1)), deparse(substitute(dmrs2))) }
  cols = match(colNames, names(dmrs1));
  par(mfrow=mfrow);
  for (i in cols) {
    cat(paste(names(dmrs1)[i], '... ', sep=''));
    if (abs) {
      d1 = abs(dmrs1[,i]); d2 = abs(dmrs2[,i]);
    } else {
      d1 = dmrs1[,i]; d2 = dmrs2[,i];
    }
    pval = signif(wilcox.test(d1, d2)$p.value, 3);
    boxplot(d1, d2, main=paste(names(dmrs1)[i], ': p=',pval, sep=''), frame.plot=F, names=grpNames, ...);
  }
}

################################################################################################
#### load A.burtoni gene annotation information from multiple files and organize it
## Args
##  dir: name of directory holding all files
## Output
##  list with 4 elements:
##   1- dataframe representing a reduced version of the burtoni gff3 file from ncbi
##       also has columns with gene lengths, %GC, number of isoforms, and transcript lengths
##   2- GRanges version of 1
##   3- dataframe with results from BLASTing all burtoni transcripts against other fish 
##       includes human homolog entrez and ensembl ids and GO terms, all from BioMart
##   4- vector of LOCs for burtoni genes the Fernald lab has a priori interest in
## Comment
##  should add info here about scripts used to generate all the files

.loadBurtoniGenes = function (dir) {
  annoDIR = dir;
  gffFILE = 'ref_AstBur1.0_scaffolds.clean.translate.final.combo.gff3_bedtools_nucLOC';
  gffexonFILE = 'ref_AstBur1.0_scaffolds.clean.translate.final.gff3_exonLengths'
  recipBlastFILE = 'H_burtoni_rna_blastx_FISH_ENS_top1_reciprocalBackFrom_Drer_Olat_Onil_Trub_ENS_pep_noComments_passRecipHsENS_withGO';
  goisFILE = 'Fernald.DAVID - Sheet1.csv';
  
  # gff files with gene LOCs, length, %GC, num isoforms, avg and max transcript length
  gff = read.table(paste(annoDIR, '/', gffFILE, sep=''), sep='\t', header=F);
  gff_exonLength0 = read.table(paste(annoDIR, '/', gffexonFILE, sep=''), sep='\t', header=T, row.names=1);
  
  # results of reciprocal blast to other fish, with hs homolog info and GO terms, "ANNO"
  anno = get(load(paste(annoDIR, '/', recipBlastFILE, sep=''))); rm(ANNO)   # loads object named ANNO
  
  # Fernald lab genes of interest from googleDoc
  gois0 = read.csv(paste(annoDIR, '/', goisFILE, sep=''));
  gois = paste('LOC',as.vector(na.omit(as.vector(as.matrix(gois0[, grepl('LOC',names(gois0))])))),sep='');
  
  gffSplit = unlist(strsplit(gff$V9, ';'));
  gffLOCs = grep('gene=', gffSplit);
  gffLOCs = unlist(strsplit(gffSplit[gffLOCs], 'gene='));   # use gsub
  gffLOCs = gffLOCs[seq(2,length(gffLOCs),2)];
  rownames(gff) = gffLOCs;
  gff2 = cbind(gff, gff_exonLength0[match(rownames(gff), rownames(gff_exonLength0)) ,]);
  names(gff2)[10:14]=c('%GC','len.gene','num.isoforms','avg.len.transcript','max.len.transcript');
  gff2gr = GRanges(seqnames=gff2$V1, ranges=IRanges(start=gff2$V4,end=gff2$V5), strand=gff2$V7, mcols=rownames(gff2));
  
  return(list(gff=gff2, gffGR=gff2gr, recipBL=anno, gois=gois));
  
}

.loadBurtoniGenesApr2016 = function () {
  load('/Users/abseq/Documents/_annotationsDec2015/WORKSPACE_parseGffFileToLookup_betterScaffoldTranslations_Feb2016.RData');
  
  # add biotype column to gffGenes and translate NW_ ids to scaffold_ ids
  geneSyms = .parseGffMetaCol(gffGenes, pattern='gene=');
  geneIDs = .parseGffMetaCol(gffGenes, pattern='Dbxref=GeneID:');
  geneBiotypes = .parseGffMetaCol(gffGenes, pattern='gene_biotype=');
  gffGenes = as.data.frame(cbind(gffGenes, biotype=geneBiotypes));
  
  # get full gene names
  #gffProduct = subset(gff, grepl('product=', V9));
  #geneSymsProduct = .parseGffMetaCol(gffProduct, pattern='gene=');
  #geneProducts = .parseGffMetaCol(gffProduct, pattern='product=');
  
  # need to translate from RefSeq scaffold ids (NW_) to scaffold numbers from BROAD (scaffold_)
  # set comment.char='' since the header line starts with a #
  scaffoldTable = read.table('/Users/abseq/Documents/_annotationsDec2015/scaffold_names', header=T, sep='\t', comment.char='', stringsAsFactors=F);
  # need to change format from, e.g. scaffold00001 to scaffold_1
  # also need to subtract one to make scaffold nums 0-based to match BROAD
  scaffoldNums = as.vector(na.omit(as.numeric(unlist(strsplit(scaffoldTable$Genome.Center.name, 'scaffold')))));
  scaffoldNums = paste0('scaffold_', scaffoldNums-1);
  scaffoldMap0 = as.data.frame(cbind(RefSeq=scaffoldTable$RefSeq.Accession.version, BROAD=scaffoldNums));
  
  scaffoldMap = scaffoldMap0[match(gffGenes$V1, scaffoldMap0$RefSeq), ];
  if (all(scaffoldMap$RefSeq==gffGenes$V1)) { gffGenes$V1 = scaffoldMap$BROAD }
  
  scaffoldMap2 = scaffoldMap0[match(gff$V1, scaffoldMap0$RefSeq), ];
  if (all(scaffoldMap2$RefSeq==gff$V1)) { gff$V1 = scaffoldMap2$BROAD }
  
  # make genomics ranges object out of gffGenes
  mcols = data.frame(geneSym=geneSyms, geneID=geneIDs, biotype=gffGenes$biotype, GC=gffGenes$V11, isoforms=gffGenes$numTranscripts);
  gffGenes2gr = GRanges(seqnames=gffGenes$V1, 
                        ranges=IRanges(start=gffGenes$V4,end=gffGenes$V5), 
                        strand=gffGenes$V7, 
                        mcols=mcols);
  
  # get typical Fernald lab genes of interest
  gois0 = read.csv('/Volumes/fishstudies-1/_Burtoni_annotations/Fernald.DAVID - Sheet1.csv');
  goiIDs = as.vector(na.omit(as.vector(as.matrix(gois0[, grepl('LOC',names(gois0))]))));
  
  # get mappings of burtoni genes to human homologs
  anno = read.table('/Users/abseq/Documents/_annotationsDec2015/passedReciprocalBlastp_wHomologs_v2_plusHomologsForFailedWithSymbols.tsv', sep='\t', header=T, quote='')
  
  gffGenesDF = as.data.frame(gffGenes2gr);
  names(gffGenesDF) = gsub('mcols.','',names(gffGenesDF),fixed=T);
  
  # put gffGenes dataframe and GR object into list with the reciprocal blastp lookup table (with homologs) and gois
  an = list(gffGenesDF=gffGenesDF, gffGenesGR=gffGenes2gr, gff=gff, recipBL=anno, gois=goiIDs, lookup=lookup);
  
  return(an);
}

################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
##  

# .makeGRForGffFeature = function (gff, fCol, startCol, endCol) {
#   
# }

################################################################################################
#### 
## Args
##  dir: name of directory holding all files
## Output
##  list with 2 elements:
##   1- list with 3 dataframes
##   2- list with 3 GRanges objects
## Comment
##  should add info here about scripts used to generate all the files

.loadBurtoniOtherFeatures = function (dir) {
  annoDIR = dir;
  lncFILE = 'abur.lnc.final.gtf';
  teFILE = 'Abur_final_TE.bed';
  miFILE = 'abur_miRNAs-130326.fix.bed';
  snpFILE = 'Assembly_SNPs.noHeader.gff3';
  df = list(lnc=read.table(paste(annoDIR, '/', lncFILE, sep=''), sep='\t', header=F,stringsAsFactors=F),
            te=read.table(paste(annoDIR, '/', teFILE, sep=''), sep='\t', header=F,stringsAsFactors=F),
            mi=read.table(paste(annoDIR, '/', miFILE, sep=''), sep='\t', header=F,stringsAsFactors=F),
            snp=read.table(paste(annoDIR, '/', snpFILE, sep=''), sep='\t', header=F,stringsAsFactors=F)
            );
  cat(paste('lncRNAs... ',sep=''));
  tmp = unlist(strsplit(df$lnc$V9, '; '));
  tmp = gsub('gene_id ', '', tmp[grepl('gene',tmp)]);
  lncGR = GRanges(seqnames=df$lnc$V1, ranges=IRanges(start=df$lnc$V4, end=df$lnc$V5), strand=df$lnc$V7, mcols=tmp);
  cat(paste('TEs... ',sep=''));
  tmp = apply(df$te, 1, function(f) gsub(' ','',paste(f[4], '_', f[1], ':', f[2], '-', f[3], sep='')));
  teGR = GRanges(seqnames=df$te$V1, ranges=IRanges(start=df$te$V2, end=df$te$V3), strand=df$te$V6, mcols=tmp);
  cat(paste('miRNAs... ',sep=''));
  miGR = GRanges(seqnames=df$mi$V1, ranges=IRanges(start=df$mi$V2, end=df$mi$V3), strand=df$mi$V6, mcols=df$mi$V4);
  cat(paste('SNPs... ',sep=''));
  snpGR = GRanges(seqnames=df$snp$V1, ranges=IRanges(start=df$snp$V4, end=df$snp$V5), strand='*', mcols=df$snp$V9);
  gr = list(lnc=lncGR, te=teGR, mi=miGR, snp=snpGR);
  return(list(df=df, gr=gr));
}

################################################################################################
#### annotate DMRs from dmrFinder() with A.burtoni annotation info in new columns
## Args
##  DMRs: dataframe of DMRs from dmrFinder()
##  ABgenesGR: GRanges version of burtoni gene gff, $gffGR from .loadBurtoniGenes()
##  maxDist: arg to findOverlaps(), max distance between DMR and gene to be considered overlapping 
##  strand: right now hardcoded to ignoer strand, need to add functionality
##  withHomologsAndGO: logical indicating whether to also add human homolog info and GO terms
##  recipBL: $recipBL from .loadBurtoniGenes(), only needed if withHomologsAndGO=TRUE
## Output
##  dataframe of input DMRs with additional columns holding annotation information
## Comment
##  this only works specifically with the annotation files loaded by .loadBurtoniGenes()

# ABgenesGR = should be $gffGR from .loadBurtoniGenes()
# .addGenesToDMRs = function (DMRs, ABgenesGR, maxDist=2000, strand=NULL, withHomologs=T, recipBL=NULL, recipCols=NULL, ...) {
#   DMRs.gr = GRanges(seqnames=DMRs$chr, ranges=IRanges(start=DMRs$start, end=DMRs$end), strand='*'); # strand
#   OV = as.data.frame(findOverlaps(DMRs.gr, ABgenesGR, ignore.strand=T, maxgap=maxDist, ...));          # strand
#   fillNA = rep(NA,nrow(DMRs));
#   DMRsOV = cbind(DMRs, gene=fillNA, chr.gene=fillNA, start.gene=fillNA, end.gene=fillNA, strand.gene=fillNA);
#   OVspl = split(OV, OV$queryHits);
#   for (i in 1:length(OVspl)) {
#     drow = OVspl[[i]]$queryHits[1];
#     grows = OVspl[[i]]$subjectHits;
#     hits = mcols(ABgenesGR)[grows, 1];
#     if ((length(hits) > 1) & (length(unique(hits)) == 1)) { hits = unique(hits) }
#     DMRsOV$gene[drow] = paste(hits, collapse=',');     # will be issues parsing gene descriptions if use ','
#     if ((length(hits) == 1) & (length(grows) == 1)) {
#       tmp = ABgenesGR[grows]
#       DMRsOV$chr.gene[drow] = as.character(seqnames(tmp));
#       DMRsOV$start.gene[drow] = data.frame(ranges(tmp))$start;
#       DMRsOV$end.gene[drow] = data.frame(ranges(tmp))$end;
#       DMRsOV$strand.gene[drow] = as.character(strand(tmp));
#     } else {
#       next;
#     }
#   }
#   if (withHomologs & is.data.frame(recipBL) & is.numeric(recipCols)) {
#     return(.addHsHomologsToDMRsWithGenes(DMRsOV, recipBL, recipCols)$dmrs);
#   }
#   return(DMRsOV);
# }

.addGenesToDMRsNew = function (DMRs, ABgenesGR, maxDist=2000, strand=NULL, withHomologs=F, recipBL=NULL, recipCols=NULL, ...) {
  DMRs.gr = GRanges(seqnames=DMRs$chr, ranges=IRanges(start=DMRs$start, end=DMRs$end), strand='*'); # strand
  OV = as.data.frame(findOverlaps(DMRs.gr, ABgenesGR, ignore.strand=T, maxgap=maxDist, ...));          # strand
  fillNA = rep(NA,nrow(DMRs));
  DMRsOV = cbind(DMRs, gene=fillNA, chr.gene=fillNA, start.gene=fillNA, end.gene=fillNA, strand.gene=fillNA, type.gene=fillNA);
  OVspl = split(OV, OV$queryHits);
  for (i in 1:length(OVspl)) {
    drow = OVspl[[i]]$queryHits[1];
    grows = OVspl[[i]]$subjectHits;
    hits = mcols(ABgenesGR)[grows, 1];
    if ((length(hits) > 1) & (length(unique(hits)) == 1)) { hits = unique(hits) }
    DMRsOV$gene[drow] = paste(hits, collapse=',');     # will be issues parsing gene descriptions if use ','
    if ((length(hits) == 1) & (length(grows) == 1)) {
      tmp = ABgenesGR[grows]
      DMRsOV$chr.gene[drow] = as.character(seqnames(tmp));
      DMRsOV$start.gene[drow] = data.frame(ranges(tmp))$start;
      DMRsOV$end.gene[drow] = data.frame(ranges(tmp))$end;
      DMRsOV$strand.gene[drow] = as.character(strand(tmp));
      DMRsOV$type.gene[drow] = as.character(mcols(tmp)[,3]);     #  careful of column
    } else {
      next;
    }
  }
  if (withHomologs & is.data.frame(recipBL) & is.numeric(recipCols)) {
    return(.addHsHomologsToDMRsWithGenes(DMRsOV, recipBL, recipCols)$dmrs);
  }
  return(DMRsOV);
}

################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
##  functionality available by setting withHomologsAndGO=T when running .addGenesToDMRs()

.addHsHomologsToDMRsWithGenes = function (DMRsOV, recipBL, recipCols) {
  multiRows = which(sapply(strsplit(DMRsOV$gene, ','), length) > 1);
  tmpAnno = recipBL[match(DMRsOV$gene, recipBL$gene), ];
  tmpAnno = cbind(DMRsOV, tmpAnno[,recipCols]);
  cols = seq(from=(ncol(tmpAnno)-length(recipCols)+1), to=ncol(tmpAnno), by=1); 
  for (i in multiRows) {
    genes = unlist(strsplit(tmpAnno$gene[i], ','));
    hits = recipBL[match(genes, recipBL$gene), recipCols];
    tmpAnno[i, cols] = apply(hits, 2, function(f) paste(f, collapse=','));
  }
  return(list(dmrs=tmpAnno, multiRows=multiRows));
}

################################################################################################
#### 
## Args
##  
## Output
##  
## Comment

.quickGO = function (gois, recipBL, thresh=.05) {
  goiTerms = recipBL[match(gois, an$recipBL$gene), ]$go_id;
  goiTerms = as.vector(na.omit(unlist(strsplit(goiTerms, '\t'))));
  allTerms = as.vector(na.omit(unlist(strsplit(recipBL$go_id, '\t'))));
  termnames = names(table(goiTerms));
  for (t in 1:length(termnames)) {
    ygoi = sum(goiTerms==termnames[t]);
    yall = sum(allTerms==termnames[t]);
    mat = matrix(c(ygoi, length(gois)-ygoi, yall, nrow(recipBL)-yall), ncol=2);
    res = fisher.test(mat);
    if (res$p.value < thresh & ygoi>2) {
      print(termnames[t]);
      print(mat)
      print(res)
    }
  }
}

################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
##  this only works specifically with the annotation files loaded by .loadBurtoniOtherFeatures()

.addOtherFeaturesToDMRs = function (DMRs, ABfeaturesGR, colname=NULL, maxDist=2000, minOv=NULL, strand=NULL, ...) {
  if (is.null(minOv)) { minOv = 1 }
 # DMRlengths = abs(DMRs$start - DMRs$end);
  if (maxDist > 0) {minOv = 0}
  DMRs.gr = GRanges(seqnames=DMRs$chr, ranges=IRanges(start=DMRs$start, end=DMRs$end), strand='*'); # strand
  OV = as.data.frame(findOverlaps(DMRs.gr, ABfeaturesGR, ignore.strand=T, maxgap=maxDist, minoverlap=minOv, ...));          # strand
  fillNA = rep(NA,nrow(DMRs));
  colname = ifelse(test=is.null(colname), 
                   yes=deparse(substitute(ABfeaturesGR)),
                   no=as.character(colname)
                   );
  DMRsOV = cbind(DMRs, fillNA, fillNA, fillNA, fillNA, fillNA);
  names(DMRsOV)[ncol(DMRs)+1] = colname;
  names(DMRsOV)[(ncol(DMRs)+2):(ncol(DMRsOV))] = paste(c('chr','start','end','strand'), colname, sep='.');
  OVspl = split(OV, OV$queryHits);
  for (i in 1:length(OVspl)) {
    drow = OVspl[[i]]$queryHits[1];
    grows = OVspl[[i]]$subjectHits;
    hits = mcols(ABfeaturesGR)[grows, 1];
    if ((length(hits) > 1) & (length(unique(hits)) == 1)) { hits = unique(hits) }
    DMRsOV[drow, ncol(DMRs)+1] = paste(hits, collapse=',');
    if ((length(hits) == 1) & (length(grows) == 1)) {
      tmp = ABfeaturesGR[grows]
      chr = as.character(seqnames(tmp));
      start = data.frame(ranges(tmp))$start;
      end = data.frame(ranges(tmp))$end;
      strand = as.character(strand(tmp));
      DMRsOV[drow, (ncol(DMRs)+2):(ncol(DMRsOV))] = c(chr, start, end, strand);
    } else {
      next;
    }
  }
  return(list(dmrs=DMRsOV,ov=OV));
}

################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
##  

.dmrDistFromFeatureStart = function (dmrs, startdmr=2, enddmr=3, startfeature=19, endfeature=20, strandfeature=21) {
  dmrs = as.data.frame(cbind(dmrs, dist=rep(NA, nrow(dmrs))));
  for (i in 1:nrow(dmrs)) {
    if (is.na(dmrs[i,strandfeature])) {
      next;
    } else if (dmrs[i,strandfeature] == '+') {
      dmrs$dist[i] = (dmrs[i, startfeature] - dmrs[i, startdmr]);
    } else if (dmrs[i,strandfeature] == '-') {
      dmrs$dist[i] = (dmrs[i, endfeature] - dmrs[i, enddmr]);
    } else {
      stop(paste0('stopped on row ', i));
    }
  }
  return(dmrs);
}



################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
##  

.getMultiRows = function (df, colname='gene', sep=',', num=1) {
  col = match(colname, names(df));
  multiRows = which(sapply(strsplit(df[,col], split=sep), length) > num);
  return(multiRows);
}

################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
## 

.buildGeneBackgroundForGO = function (recipBL, idType='entrezgene', allowMulti=2, multiDelim='\t', write=T, prefix=NULL) {
  idCol = match(idType, names(recipBL));
  genes00 = recipBL[, idCol];
  genes0 = genes00[(!is.na(genes00)) & (genes00 != 'NA')];
  multi = grepl(multiDelim, genes0);
  genes1 = genes0[!multi];
  keep = c();
  for (i in which(multi)) {
    spl = strsplit(genes0[i], multiDelim, fixed=T)[[1]]; 
    if (length(spl) <= allowMulti) {
      keep = c(keep, spl)
    } else {
      keep = c(keep, spl[1:allowMulti])
    }
  }
  genes = unique(c(genes1, keep));
  genes = genes[genes != 'NA'];
  if (write) {
    file = paste('BG_', idType, '_multi', allowMulti, '.txt', sep='');
    if (is.character(prefix)) { file = paste(prefix, file, sep='_') }
    write.table(genes, file=file, row.names=F, col.names=F, quote=F, sep='\t');
  }
  return(genes);
}


################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
##  could be faster if vectorized but loop is easier for now

.getIDsForGOFromDMRs = function (DMRs, ABgeneColName='gene', idType='entrezgene', allowMulti=2, write=F, prefix=NULL) {
  idCol = match(idType, names(DMRs));
  abCol = match(ABgeneColName, names(DMRs));
  genes = c();
  for (i in 1:nrow(DMRs)) {
    this_id = DMRs[i, idCol];
    ab_id = DMRs[i, abCol];
    if (is.na(this_id)) { next }
    multicheck = grepl('\t|,', this_id);
    if (multicheck) {
      spl0 = unlist(strsplit(this_id, '\t|,'));
      abspl = unlist(strsplit(ab_id, '\t|,'));
      spl = unique(spl0[!is.na(spl0) & (spl0 != 'NA')]);
      thresh = allowMulti * length(abspl);
      if (length(spl) <= thresh) {
        genes = c(genes, spl);
      } else {
        genes = c(genes, spl[1:thresh]);
      }
    } else {
      genes = c(genes, this_id);
    }
  }
  genes = unique(gsub(' ', '', genes));
  if (write) {
    file = paste(idType, '_multi', allowMulti, '.txt', sep='');
    if (is.character(prefix)) { file = paste(prefix, file, sep='_') }
    write.table(genes, file=file, row.names=F, col.names=F, quote=F, sep='\t');
  }
  return(genes);
}

################################################################################################
#### get IDs for genes near DMRs
## Args
##  DMRsOV:
##  
##  IDtype:
##  allowMulti:
##  write:
##  
## Output
##  if write=T, writes text files to working directory
##  also list with 3 elements:
##   1-
##   2-
##   3-
## Comment
##

.getIDsForGOFromDMRsAcrossGroups = function (DMRsOV, ABgeneColName='gene', idType='entrezgene', allowMulti=2, write=T, prefix=NULL) {
  # build lists for hyper and hypo DMRs
  idColDMR = match(idType, names(DMRsOV));
  for (GRP in c('hypo','hyper')) {
    DMRs = subset(DMRsOV, direction==GRP);
    fprefix = ifelse(is.character(prefix), paste(prefix, GRP, sep='_'), GRP);
    assign(GRP, .getIDsForGOFromDMRs(DMRs, ABgeneColName=ABgeneColName, idType=idType, allowMulti=allowMulti, write=write, prefix=fprefix));
  }
  allIDs = c(hypo, hyper)
  if (write) {
    fileBase = paste(idType, '_multi', allowMulti, '.txt', sep='');
    for (f in c('hypo', 'hyper', 'allIDs')) {
      thisfile = paste('file', f, sep='');
      if (is.character(prefix)) {
        assign(thisfile, paste(prefix, f, fileBase, sep='_'));
      } else {
        assign(thisfile, paste(f, fileBase, sep='_'));
      }
      write.table(get(f), file=get(thisfile), row.names=F, col.names=F, quote=F, sep='\t')
    }
  }
  return(list(hypo=hypo, hyper=hyper, all=allIDs));
}
 
################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
##

.GOFunctionAllOntologies = function (genes, refGenes, ...) {
  cat('\n==========  BP  ==========\n');
  tmpBP = GOFunction(interestGenes=genes, refGenes=refGenes, ontology='BP', ...);print(tmpBP)
  if (!is.null(tmpBP)) {
    tmpBP = as.data.frame(cbind(ontology=rep('BP',nrow(tmpBP)), tmpBP));
  }
  cat('\n==========  CC  ==========\n');
  tmpCC = GOFunction(interestGenes=genes, refGenes=refGenes, ontology='CC', ...);print(tmpCC)
  if (!is.null(tmpCC)) {
    tmpCC = as.data.frame(cbind(ontology=rep('CC',nrow(tmpCC)), tmpCC));
  }
  cat('\n==========  MF  ==========\n');
  tmpMF = GOFunction(interestGenes=genes, refGenes=refGenes, ontology='MF', ...);print(tmpMF)
  if (!is.null(tmpMF)) {
    tmpMF = as.data.frame(cbind(ontology=rep('MF',nrow(tmpMF)), tmpMF));
  }
  tmpAll = as.data.frame(matrix(ncol=7, dimnames=list(1, c('ontology', 'goid', 'name', 'refnum', 'interestnum', 'pvalue', 'adjustp'))))
  for (res in c('tmpBP','tmpCC','tmpMF')) {
    if (!is.null(get(res))) {
      tmpAll = as.data.frame(rbind(tmpAll, get(res)))
    }
  }
  tmpAll = tmpAll[-1, ];
  return(tmpAll[order(tmpAll$adjustp, -tmpAll$interestnum), ]);
}

.enrichGOAllOntologies = function (genes, refGenes, ...) {
  cat('\n==========  BP  ==========\n');
  tmpBP = summary(enrichGO(genes, universe=refGenes, ont='BP', ...));
  if (!is.null(tmpBP)) {
    tmpBP = as.data.frame(cbind(Ontology=rep('BP',nrow(tmpBP)), tmpBP));print(dim(tmpBP))
  }
  cat('\n==========  CC  ==========\n');
  tmpCC = summary(enrichGO(genes, universe=refGenes, ont='CC',  ...));
  if (!is.null(tmpCC)) {
    tmpCC = as.data.frame(cbind(Ontology=rep('CC',nrow(tmpCC)), tmpCC));print(dim(tmpCC))
  }
  cat('\n==========  MF  ==========\n');
  tmpMF = summary(enrichGO(genes, universe=refGenes, ont='MF', ...));
  if (!is.null(tmpMF)) {
    tmpMF = as.data.frame(cbind(Ontology=rep('MF',nrow(tmpMF)), tmpMF));print(dim(tmpMF))
  }
  tmpAll = as.data.frame(matrix(ncol=10, dimnames=list(1, c('Ontology', 'ID', 'Description', 'GeneRatio', 'BgRatio', 'pvalue', 'p.adjust', 'qvalue','geneID','Count'))))
  for (res in c('tmpBP','tmpCC','tmpMF')) {
    if (!is.null(get(res))) {
      this = get(res); 
      tmpAll = as.data.frame(rbind(tmpAll, this))
    }
  }
  tmpAll = tmpAll[-1, ];
  return(tmpAll[order(tmpAll$qvalue), ]);
}



################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
##

.geneGffToGRanges = function (gff, seqnamesCol, startCol, endCol, strandCol, mcols='rownames') {
  if (!all(is.numeric(c(seqnamesCol, startCol, endCol, strandCol)))) { stop('*Col args must be numeric') }
  if (mcols == 'rownames') {
    mcols = rownames(gff);
  } else if (is.numeric(mcols)) {
    mcols = gff[, mcols];
  } else if (!is.null(mcols)) {
    
  }
  gr = GRanges(seqnames=gff[,seqnamesCol], 
               ranges=IRanges(start=gff[,],end=gff[,]), 
               strand=gff[,], mcols=rownames(gff)
               );
}


################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
##

.parseGffMetaCol = function(gff, metaCol=9, delim=';', pattern='Dbxref=GeneID:') {
  metaSplit = strsplit(gff[, metaCol], delim);
  check = sapply(metaSplit, 
                 function(f) any(grepl(pattern, f)));
  if (!all(check) | sum(check)!=nrow(gff)) {
    stop('all lines in gff do not contain arg pattern') 
  }
  res = grep(pattern, unlist(metaSplit), value=TRUE);
  return(gsub(pattern, '', res));
}

.buildFeatureGR = function (gff, feature, chrCol=1, startCol=4, endCol=5, featureCol=3, idPattern='gene=', metaCol=9, delim=';', strandCol=7) {
  cat('Getting gene ids... ');
  ids = .parseGffMetaCol(gff=gff, metaCol=metaCol, delim=delim, pattern=idPattern);
  fRows = gff[, featureCol] == feature;
  feat0 = paste(gff[fRows, chrCol], gff[fRows, startCol], gff[fRows, endCol], gff[fRows, strandCol]);
  names(feat0) = ids[fRows];
  feat = feat0[!duplicated(feat0)];
  cat('Building GRanges... seqnames...');
  fseqnames = unlist(lapply(strsplit(feat, ' '), 
                            function(f) f[1]));
  cat('starts...');
  fstarts = as.numeric(unlist(lapply(strsplit(feat, ' '), 
                                     function(f) f[2])));
  cat('ends...');
  fends = as.numeric(unlist(lapply(strsplit(feat, ' '), 
                                   function(f) f[3])));
  cat('strand...');
  fstrand = unlist(lapply(strsplit(feat, ' '), 
                          function(f) f[4]));
  featGR = GRanges(seqnames=fseqnames, 
                   ranges=IRanges(start=fstarts, end=fends), 
                   mcols=names(feat), 
                   strand=fstrand);
  return(featGR);
}




################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
##

.getLASAGNAResults = function (tsvfile, dmrs) {
  tfs = read.table(tsvfile, header=F, sep='\t', fill=T, quote='');
  scaffoldRows = grep('^scaffold', tfs$V1);
  tflist = list()
  for (i in 1:length(scaffoldRows)) {
    if (i == length(scaffoldRows)) {
      this_dmr = tfs[(scaffoldRows[i]+2):nrow(tfs), ];
    } else {
      this_dmr = tfs[(scaffoldRows[i]+2):(scaffoldRows[i+1]-1), ]; 
    }
    this_dmrName = tfs[(scaffoldRows[i]), 1];
    this_chr = unlist(strsplit(this_dmrName, ':'))[1];
    this_st_end = as.numeric(unlist(strsplit(unlist(strsplit(this_dmrName, ':'))[2], '-')));
    this_strand = dmrs[dmrs$chr==this_chr & dmrs$start==this_st_end[1] & dmrs$end==this_st_end[2], ]$strand.gene;
    this_dmr = subset(this_dmr, V4==this_strand)
    tflist[[i]] = this_dmr;
    names(tflist)[i] = this_dmrName
  }
  return(tflist);
}

################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
##

.getExonStarts = function (gfffile, featureCol=3, startCol=4, endCol=5, strandCol=7, frameCol=8) {
  cat(paste0('Reading ', gfffile, '... \n'));
  gff = read.table(gfffile,header=F,sep='\t');
  gff = gff[grepl('gene=LOC', gff$V9), ];
  gffSplit = unlist(strsplit(gff$V9, ';'));
  gffLOCs = grep('gene=LOC', gffSplit);
  gffLOCs = gsub('gene=', '', gffSplit[gffLOCs]);
  cat(paste0('Splitting gff by LOC ids... \n'));
  gffSplit = split(gff, gffLOCs);
  cat(paste0('Finding exon/CDS starts... \n'));
  res = as.data.frame(matrix(nrow=length(gffSplit), ncol=4, dimnames=list(names(gffSplit),c('chr','utr','cds','strand'))));
  for (i in 1:length(gffSplit)) {
    thisLOC = gffSplit[[i]];
    strand = thisLOC[1, strandCol];
    exonRows = thisLOC[,featureCol]=='exon'
    exonPos = c(thisLOC[exonRows, startCol], thisLOC[exonRows, endCol]);
    CDSRows = thisLOC[,featureCol]=='CDS';
    if (sum(CDSRows) == 0) {
      if (strand == '+') {
        res$utr[i] = min(exonPos);
        res$cds[i] = NA;
        res$strand[i] = strand;
        res$chr[i] = thisLOC$V1[1];
      } else if (strand == '-') {
        res$utr[i] = max(exonPos);
        res$cds[i] = NA;
        res$strand[i] = strand;
        res$chr[i] = thisLOC$V1[1];
      } else {
        stop();
      }
    } else {
      CDSPos = c(thisLOC[CDSRows, startCol], thisLOC[CDSRows, endCol]);
      if (strand == '+') { 
        utr_start = min(exonPos);
        cds_start = min(CDSPos);
        frame = suppressWarnings(unique(as.numeric(thisLOC[thisLOC[,startCol]==cds_start, frameCol]))); 
        frame = frame[!is.na(frame)];
        if (length(unique(na.omit(frame))) > 1){stop()}
        cds_start = cds_start + frame;
      } else if (strand == '-') { 
        utr_start = max(exonPos);
        cds_start = max(CDSPos);
        frame = suppressWarnings(unique(as.numeric(thisLOC[which(thisLOC[,endCol]==cds_start), frameCol]))); 
        frame = frame[!is.na(frame)];
        if (length(unique(na.omit(frame))) > 1){stop()}
        cds_start = cds_start - frame;
      } else {
        stop()
      }
      res$utr[i] = utr_start;
      res$cds[i] = cds_start;
      res$strand[i] = strand;
      res$chr[i] = thisLOC$V1[1];
    }
  }
  res$utr[is.infinite(res$utr)] = NA;
  cat(paste0('Building output... \n'));
  res2 = res[!is.na(res$utr) & !is.na(res$cds), ];
  res2forGR = res2;
  names(res2forGR)[2:3] = c('start','end');
  for (i in 1:nrow(res2forGR)) {
    if (res2forGR$strand[i] == '-') {
      st = res2forGR$end[i];
      en = res2forGR$start[i];
      res2forGR$start[i] = st;
      res2forGR$end[i] = en;
    } 
  }
  res2GR = GRanges(seqnames=res2forGR$chr, ranges=IRanges(start=res2forGR$start, end=res2forGR$end), strand=res2forGR$strand);
  mcols(res2GR) = rownames(res2forGR)
  
  return(list(df=res, dfForGR=res2forGR, gr=res2GR));
}

################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
##

#org.Hs.egGO2ALLEGSlist = as.list(org.Hs.egGO2ALLEGS);
.getGenesForGOTerms = function (terms, maps=org.Hs.egGO2ALLEGSlist, genes=NULL) {
  glist = list();
  for (t in terms) {
    g = maps[[t]];
    if (!is.null(genes)) {
      glist[[t]] = intersect(genes, g);
    } else {
      glist[[t]] = g;					### may include duplicates
    }
  }
  return(glist);
}

# x3 = .getGenesForGOTerms(terms=gofCC$turquoise$goid, genes=modGenes.hsEntrez$turquoise[,2])
# x3ENS = lapply(x3, function(f) modGenes.hsEntrez$turquoise[modGenes.hsEntrez$turquoise[,2] %in% f, 1]);

# gene ids should be Entrez
.addGenesToGOFunctionResults = function (GOFunctionAllOntologiesOutput, maps=org.Hs.egGO2ALLEGSlist, modulegenes) {
  genelist = .getGenesForGOTerms(GOFunctionAllOntologiesOutput$goid, genes=modulegenes);
  new = as.data.frame(cbind(GOFunctionAllOntologiesOutput, genes=rep(NA, nrow(GOFunctionAllOntologiesOutput))));
  for (row in 1:nrow(new)) {
    tgenes = genelist[[match(new$goid[row], names(genelist))]];
    new$genes[row] = paste(tgenes, collapse=',');
  }
  return(new);
}

#tmp = .addGenesToGOFunctionResultsList(tmp, modGenesList=lapply(modGenes.hsEntrez, function(f) f[,2]))

.addGenesToGOFunctionResultsList = function (GOFunctionModulesOutput, maps=org.Hs.egGO2ALLEGSlist, modGenesList) {
  newlist = GOFunctionModulesOutput;
  for (m in 1:length(newlist)) {
    if (nrow(newlist[[m]]) == 0) {
      cat('skipping ', names(newlist)[m], ' since no terms\n', sep='');
      next;
    }
    newlist[[m]] = .addGenesToGOFunctionResults(newlist[[m]], modulegenes=modGenesList[[m]]);
  }
  return(newlist);
}

.addNumToGOFunctionResults = function (GOFunctionAllOntologiesOutput, nums, entrezToENS=NULL, f='median', colname=NULL) {
  if (any(names(GOFunctionAllOntologiesOutput)=='genes')) {
    new = as.data.frame(cbind(GOFunctionAllOntologiesOutput, rep(NA, nrow(GOFunctionAllOntologiesOutput))));
    for (row in 1:nrow(new)) {
      genes = unlist(strsplit(new$genes[row], ','));
      if (!is.null(entrezToENS)) {
        genes = entrezToENS[match(genes, entrezToENS[,2]), 1];
      }
      new[row, ncol(new)] = .getSomeNumForGenes(genes=genes, nums=nums, f=f);
    }
  } else {
    stop('add genes to results first');
  }
  if (is.character(colname)) {
    names(new)[ncol(new)] = colname;
  } else {
    names(new)[ncol(new)] = paste(deparse(substitute(nums)), '.', f, sep='');
  }
  return(new);
}

# maplist should be list of vectors
# names should be burtoni gene names and elements are vectors of human entrez ids
.mapAbGenesHsEntrez = function (ids, maplist, inputType=c('hs','ab')) {
  inputType = inputType[1];
  if (!(inputType %in% c('hs','ab'))) { 
    warning('Invalid inputType, assuming input is Hs entrez ids');
    inputType = 'hs';
  }
  if (inputType == 'hs') {
    return(maplist[sapply(maplist, function(f) any(f %in% ids))]);
  } else {
    return(maplist[names(maplist) %in% ids]);
  }
}

.getGOFunctionTermGenes = function (GOFunctionAllOntologiesOutput, row, geneCol=which(names(GOFunctionAllOntologiesOutput)=='genes')) {
  go = GOFunctionAllOntologiesOutput;
  if (is.character(row)) { row = which(go$name==row | go$goid==row) };
  return(unique(as.vector(na.omit(unlist(strsplit(go[row, geneCol], ','))))));
}

.getNeighborsForGOFunctionTermGenes = function (GOFunctionAllOntologiesOutput, row, geneCol, maplist, checkGenes, sGR, stranded=TRUE, ...) {
  go = GOFunctionAllOntologiesOutput;
  termgenes = .getGOFunctionTermGenes(go, row=row, geneCol=geneCol);
  abtermgenes = .mapAbGenesHsEntrez(termgenes, maplist=maplist, inputType='hs');
  abtermgenes = abtermgenes[names(abtermgenes) %in% checkGenes];
  abtermgenesGR = sGR[mcols(sGR)[,1] %in% names(abtermgenes)];   # assumes gene names in mcols(sGR)[,1]   
  if (!stranded) { strand(abtermgenesGR) <- '*' }
  nn = .getNeighborsUnstranded(abtermgenesGR, sGR)
  return(list(termGR=abtermgenesGR, neighbors=nn));
}

################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
## df1 and df2 are probably an$recipBL and handannos

.combineListsOfVectorsWithSameNames = function (l1, l2) {
  if(!all(names(l1) == names(l2))) { stop('lists must have exact same length and names') }
  for (i in 1:length(l1)) {
    new = unique(c(l1[[i]], l2[[i]]));
    new = new[!(new %in% c('NA','NaN'))];
    new = new[!is.na(new)];
    l1[[i]] = sort(new);
  }
  return(l1)
}

.mergeBurtoniGeneToHumanEntrezMappings = function (df1, df2, df1.ab, df1.hs, df2.ab, df2.hs) {
  m1 = strsplit(df1[, df1.hs], ',');
  names(m1) = df1[, df1.ab];
  m2 = strsplit(df2[, df2.hs], ',');
  names(m2) = df2[, df2.ab];
  commnames = intersect(names(m1), names(m2));
  commlist = .combineListsOfVectorsWithSameNames(m1[commnames], m2[commnames]);
  only1names = setdiff(names(m1), names(m2));
  only2names = setdiff(names(m2), names(m1));
  new = c(m1[only1names], m2[only2names], commlist);
  return(new[order(names(new))]);
}


################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
##

.makeRaggedMatrix = function (listOfVecs) {
  numvecs = length(listOfVecs);
  maxveclength = max(sapply(listOfVecs,length));
  mat = matrix(rep(NA, (maxveclength*numvecs)), ncol=numvecs);
  for (v in 1:numvecs) {
    thisvec = listOfVecs[[v]];
    if (length(thisvec)==0) { thisvec = NA }
    mat[1:length(thisvec), v] = thisvec;
  }
  return(mat);
}

.barplotStackedFromVecs = function (listOfVecs, cex.text=1, labOffset=.2, ...) {
  mat = .makeRaggedMatrix(listOfVecs);
  mids = barplot(mat, beside=F, ...);
  hvec = listOfVecs;
  for (v in 1:length(hvec)) {
    vec = hvec[[v]];
    text(mids[v]+labOffset, vec[1] / 2, 
         paste(names(vec)[1],' - ',round(vec[1],2)*100,'%',sep=''), 
         cex=cex.text);
    for (i in 2:length(vec)) {
      h = sum(vec[1:i]);
      d = h - sum(vec[1:(i-1)]);
      #text(mids[v]+labOffset, h-(d/2), names(vec)[i], cex=cex.text);
      text(mids[v]+labOffset, h-(d/2), 
           paste(names(vec)[i],' - ',round(vec[i],2)*100,'%',sep=''),
           cex=cex.text);
    }
  }
}

############### maybe sketchy results due to 0vs1-based or +/- strand
.checkIfSNPMatchesAnyCpGInDMR = function (chr.dmr, start.dmr, end.dmr, snp, snpAnnos, fitlist) {
  r = as.data.frame(ranges(fitlist[[which(names(fitlist)==chr.dmr)]]));
  rdmr = subset(r, start >= start.dmr & end <= end.dmr);
  scaffoldsnps = subset(snpAnnos, V1==chr.dmr);
  snppos = scaffoldsnps$V4[scaffoldsnps$V9 == snp];
  output = ifelse(any(rdmr$start == snppos), TRUE, FALSE);
  return(output);
}

.countSNPsOnCpGsInDMRs = function (dmrs, snpAnnos, fitlist) {
  nums = rep(NA, nrow(dmrs));
  for (i in 1:nrow(dmrs)) {
    if (!is.na(dmrs$snp[i])) {
      snpnum = length(strsplit(dmrs[i, which(names(dmrs)=='snp')], 'Name')[[1]]) - 1;
      if (snpnum==1) {
        nums[i] = as.numeric(.checkIfSNPMatchesAnyCpGInDMR(dmrs$chr[i], dmrs$start[i], dmrs$end[i], dmrs$snp[i], snpAnnos, fitlist));
      } else if (snpnum > 1) {
        tmp = c();
        these_snps = paste0('Name', unlist(strsplit(dmrs$snp[i], 'Name')));
        these_snps = gsub(',$', '', these_snps);
        for (j in 2:length(these_snps)) {
          tmp = c(tmp, .checkIfSNPMatchesAnyCpGInDMR(dmrs$chr[i], dmrs$start[i], dmrs$end[i], these_snps[j], snpAnnos, fitlist));
        }
        nums[i] = sum(tmp);
      } else {
        stop();
      }
    }
  }
  return(nums);
}
#################

################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
## 

.addUpstreamBPToGRanges = function (GRangesObj, upBP=1000) {
  tmpGR = GRangesObj;
  plusrows = which(as.vector(strand(tmpGR)=='+'));
  minusrows = which(as.vector(strand(tmpGR)=='-'));
  newranges = data.frame(ranges(tmpGR));
  newranges$start[plusrows] = newranges$start[plusrows] - upBP;
  newranges$start[plusrows][newranges$start[plusrows] < 0] = 0;
  newranges$end[minusrows] = newranges$end[minusrows] + upBP;
  newranges$width = newranges$width + upBP;
  ranges(tmpGR) = IRanges(start=newranges$start, end=newranges$end);
  
  oldranges = data.frame(ranges(GRangesObj));
  upranges = data.frame(ranges(tmpGR));
  upranges$end[plusrows] = oldranges$start[plusrows];
  upranges$start[minusrows] = oldranges$end[minusrows];
  upGR = tmpGR;
  ranges(upGR) = IRanges(start=upranges$start, end=upranges$end);
  
  return(list(new=tmpGR, upGR=upGR));
}

.extendGRangesSymmetric = function (gr, extension=1000) {
  newranges = as.data.frame(ranges(gr));
  newstarts = newranges$start - extension;
  newends = newranges$end + extension;
  newgr = GRanges(seqnames=seqnames(gr), ranges=IRanges(start=newstarts, end=newends), 
                  strand=strand(gr), mcols=mcols(gr));
  return(newgr);
}

################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
## 

.smoothAndTtestSingleScaffold = function (BSD, scaffold, NS, H, MAXGAP) {
  tmpBSD = subset(BSD, seqnames(BSD)==scaffold);
  seqlevels(tmpBSD) = scaffold;
  tmpSmooth = BSmooth(tmpBSD, ns=NS, h=H, maxGap=MAXGAP, mc.cores=4, parallelBy='sample');
  
}

################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
##  

.makeDMRsEntrezList = function (dmrs, entrezCol) {
  dmrsEntrezList = list();
  for (i in 1:nrow(dmrs)) {
    thisdmr = paste0(dmrs$chr[i],':',dmrs$start[i],'-',dmrs$end[i]);
    dmrentrez = as.vector(na.omit(unlist(strsplit(dmrs[i, entrezCol], '\t|,'))));
    dmrentrez = dmrentrez[dmrentrez != 'NA'];
    dmrsEntrezList[[thisdmr]] = dmrentrez;
  }
  return(dmrsEntrezList)
}

.getDMRsWithEntrezGenes = function (dmrs, entrezGenes, entrezCol) {
  dmrsEntrezList = .makeDMRsEntrezList(dmrs, entrezCol);
  dmrnames = paste0(dmrs$chr,':',dmrs$start,'-',dmrs$end);
  sublist = dmrsEntrezList[unlist(sapply(dmrsEntrezList, function(f) length(intersect(entrezGenes,f))>0))];
  dmrrows = dmrnames %in% names(sublist);
  return(list(dmrs=dmrs[dmrrows, ], entrezlist=sublist, dmrrows=dmrrows));
}


################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
##  tpm rownames must be LOCs

.addTPM = function (dmrs, tpm, geneCol=17, groups=c('ND','D','D','ND')) {
  dmrs = as.data.frame(cbind(dmrs, expND=rep(NA, nrow(dmrs)), expD=rep(NA, nrow(dmrs))));
  for (i in 1:nrow(dmrs)) {
    this_gene = dmrs[i, geneCol];
    if (is.na(this_gene) | grepl(',', this_gene)) {
      next;
    } else {
      exp = apply(tpm[grepl(this_gene, rownames(tpm)), ], 2, sum);
      dmrs$expND[i] = mean(exp[groups=='ND']);
      dmrs$expD[i] = mean(exp[groups=='D']);
    }
  }
  return(dmrs);
}


.addTPMup5kb = function (dmrs, tpm, geneCol=17, groups=c('ND','D','D','ND')) {
  dmrs = as.data.frame(cbind(dmrs, up5kb.expND=rep(NA, nrow(dmrs)), up5kb.expD=rep(NA, nrow(dmrs))));
  for (i in 1:nrow(dmrs)) {
    this_gene = dmrs[i, geneCol];
    if (is.na(this_gene) | grepl(',', this_gene)) {
      next;
    } else {
      exp = apply(tpm[grepl(this_gene, rownames(tpm)), ], 2, sum);
      dmrs$up5kb.expND[i] = mean(exp[groups=='ND']);
      dmrs$up5kb.expD[i] = mean(exp[groups=='D']);
    }
  }
  return(dmrs);
}


################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
##  

.computeStatsForDMRs = function (dmrs, columns) {
  dmrs = dmrs[, columns];
  out = list(statmat=NULL, numDmrs=nrow(dmrs), hypRatio=NULL);
  directionCol = which(names(dmrs) == 'direction');
  if (length(directionCol) > 0) {
    loopinds = (1:ncol(dmrs))[-directionCol];
    dircol = dmrs[, directionCol];
    out$hypRatio = sum(dircol == 'hypo') / sum(dircol == 'hyper');
  } else {
    loopinds = 1:ncol(dmrs);
  }
  dmrs = as.data.frame(apply(dmrs[,loopinds], 2, as.numeric));
  statmat = as.data.frame(matrix(nrow=6, ncol=1));
  rownames(statmat) = c('min','q1','median','mean','q3','max');
  for (i in loopinds) {
    statmat = as.data.frame(cbind(statmat, as.vector(summary(dmrs[,i]))));
    names(statmat)[ncol(statmat)] = names(dmrs)[i];
    if (any(dmrs[,i] < 0)) {
      statmat = as.data.frame(cbind(statmat, as.vector(summary(abs(dmrs[,i])))));
      names(statmat)[ncol(statmat)] = paste0(names(dmrs)[i], '.abs');
    }
  }
  out$statmat = statmat[,-1];
  return(out)
}


################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
##  

.plotDmrStatsAcrossFactor = function (dmrs, fact, ...) {
  par(mfrow=c(3,4), oma=c(0,2,5,2));
  ord = match(c('n','width','invdensity','group1.mean','group2.mean','tstat.sd','areaStat','maxStat','meanDiff'),
              names(dmrs));
  for (i in ord) {
    main=paste0('median: F=', signif(median(dmrs[!fact,i]),3),', T=', signif(median(dmrs[fact,i]), 3), '\n');
    WGCNA::verboseBoxplot(dmrs[,i], fact, 
                          ylab=names(dmrs)[i], xlab='', main=main,
                          col='grey', border='darkgrey', frame.plot=F, ...);
    if (any(dmrs[,i] < 0)) {
      main=paste0('median: F=', signif(median(abs(dmrs[!fact,i])),3),', T=', signif(median(abs(dmrs[fact,i])), 3), '\n');
      WGCNA::verboseBoxplot(abs(dmrs[,i]), fact, 
                            ylab=paste0('abs(', names(dmrs)[i] ,')'), xlab='', main=main,
                            col='grey', border='darkgrey', frame.plot=F, ...);
    }
  }
  title(deparse(substitute(fact)), outer=T);
}

################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
##  

.plotDmrStatsHistograms0 = function (dmrs, ...) {
  ord = match(c('n','width','invdensity','group1.mean','group2.mean','tstat.sd','areaStat','maxStat','meanDiff'),
              names(dmrs));
  par(mfrow=c(3,4), oma=c(2,2,2,2));
  for (i in ord) {
    hist(dmrs[,i], main=names(dmrs)[i], xlab='', col='grey', border='darkgrey', ...);
    if (any(dmrs[,i] < 0)) {
      hist(abs(dmrs[,i]), main=paste0('abs(', names(dmrs)[i] ,')'), xlab='', col='grey', border='darkgrey', ...);
    }
  }
  title(paste0(nrow(dmrs), ' DMRs'), outer=TRUE);
}

.plotDmrStatsHistogramsJpg = function (dmrs, file, res, ...) {
  jpeg(file=file, width=14, height=7, units='in', quality=100, type='quartz', res=res);
  .plotDmrStatsHistograms0(dmrs, ...)
  dev.off();
}


################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
##

.plotDmrStatsHeatmap = function (dmrs, columns, corMethod='p', col=heat.colors(200), file, res, ...) {
  jpeg(file=file, width=7, height=7, units='in', quality=100, type='quartz', res=res);
  heatmap(cor(dmrs[, columns], method=corMethod), symm=T, col=col, margins=c(7, 7), ...);
  dev.off();
}




################################################################################################
#### 
## Args
##  
## Output
##  
## Comment
## assumes that first column of mcols(sGR) holds names for the features

.getNearestNeighborsUnstrandedNoOverlap = function (qGR, sGR, verbose=100) {
  np = precede(qGR, sGR, ignore.strand=TRUE);#print(np)
  nf = follow(qGR, sGR, ignore.strand=TRUE);#print(nf)
  if (all(is.na(c(np, nf)))) {
    warning('No nearest neighbors for any range in query, returning NA');
    return(NA);
  }
  outlist = as.list(rep(NA, length(qGR)));
  bothNA = apply(cbind(np, nf), 1, function(f) all(is.na(f)))
  loopinds = which(!bothNA);#print(loopinds)
  for (n in loopinds) {#print('----------------------------------------');
    if (n %% verbose == 0) { cat(n, '...') }
    df = as.data.frame(matrix(nrow=2, ncol=7));
    dimnames(df) = list(c('left.of','right.of'), c('name','seqnames','start','end','strand','dist','qstream'));
    sinds = c(np[n], nf[n]);                                       #  print(sinds)
    oki = which(!is.na(sinds));                                     #  print(oki)
    df$name[oki] = mcols(sGR[sinds[oki]])[, 1];                     #  print(df)
    df$seqnames[oki] = as.character(seqnames(sGR[sinds[oki]]));     #  print(df)
    sGRrange = as.data.frame(ranges(sGR[sinds[oki]]));                 
    df$start[oki] = sGRrange$start;                                  #  print(df)
    df$end[oki] = sGRrange$end;                    #  print(df);    print(qGR[n]); print(sGR[sinds[oki]]); 
    df$dist[oki] = distance(qGR[n], sGR[sinds[oki]]);             # print(df)
    df$strand[oki] = as.character(strand(sGR[sinds[oki]]));
    df$qstream[1] = ifelse(df$strand[1]=='+', 'up', 'down');
    df$qstream[2] = ifelse(df$strand[2]=='-', 'up', 'down');
    outlist[[n]] = df; 
  }
  qr = as.data.frame(ranges(qGR));
  names(outlist) = paste0(as.character(seqnames(qGR)), ':', qr$start, '-', qr$end);
  return(outlist);
}

# .getNearestNeighborsStrandedNoOverlap = function (qGR, sGR, verbose=100) {
# #   outlist = as.list(rep(NA, length(qGR)));
# #   qr = as.data.frame(ranges(qGR));
# #   names(outlist) = paste0(as.character(seqnames(qGR)), ':', qr$start, '-', qr$end);
#   
#   qGRfwd = qGR[strand(qGR) == '+'];
#   qGRrev = qGR[strand(qGR) == '-'];
# 
#   upfwd = precede(qGRfwd, sGR, ignore.strand=FALSE);
#   dnfwd = follow(qGRfwd, sGR, ignore.strand=FALSE);
#   goodupfwd = !is.na(upfwd);
#   gooddnfwd = !is.na(dnfwd);
#   
#   upfwdHits = as.data.frame(cbind(as.data.frame(qGRfwd[goodupfwd]),
#                                   as.data.frame(sGR[upfwd[goodupfwd]]),
#                                   distance=distance(qGRfwd[goodupfwd],
#                                                     sGR[upfwd[goodupfwd]], 
#                                                     ignore.strand=F),
#                                   qstream=rep('up',sum(goodupfwd))));
#   dnfwdHits = as.data.frame(cbind(as.data.frame(qGRfwd[gooddnfwd]),
#                                   as.data.frame(sGR[dnfwd[gooddnfwd]]),
#                                   distance=distance(qGRfwd[gooddnfwd],
#                                                     sGR[dnfwd[gooddnfwd]], 
#                                                     ignore.strand=F),
#                                   qstream=rep('down',sum(gooddnfwd))));
#   
#   uprev = precede(qGRrev, sGR, ignore.strand=FALSE);
#   dnrev = follow(qGRrev, sGR, ignore.strand=FALSE);
#   gooduprev = !is.na(uprev);
#   gooddnrev = !is.na(dnrev);
#   
# }
# 
.getNearestNeighborsStrandedNoOverlap = function (qGR, sGR, verbose=100) {
  up = precede(qGR, sGR, ignore.strand=FALSE);
  dn = follow(qGR, sGR, ignore.strand=FALSE);
  goodup = !is.na(up);
  gooddn = !is.na(dn);
  
  upHits = as.data.frame(cbind(as.data.frame(qGR[goodup]),
                                  as.data.frame(sGR[up[goodup]]),
                                  distance=distance(qGR[goodup],
                                                    sGR[up[goodup]], 
                                                    ignore.strand=F),
                                  qstream=rep('up',sum(goodup))));
  dnHits = as.data.frame(cbind(as.data.frame(qGR[gooddn]),
                                  as.data.frame(sGR[dn[gooddn]]),
                                  distance=distance(qGR[gooddn],
                                                    sGR[dn[gooddn]], 
                                                    ignore.strand=F),
                                  qstream=rep('down',sum(gooddn))));
  
  df = as.data.frame(rbind(upHits, dnHits));
  dfsplit = split(df, paste0(df$seqnames,':',df$start,'-',df$end));
  
}

# returns list with same length as qGR
# elements for qGR ranges with no overlaps will be NA
.getOverlapsUnstranded = function (qGR, sGR, ...) {
  ov = findOverlaps(qGR, sGR, maxgap=0L, ignore.strand=TRUE, ...);
  if (length(ov) == 0) {
    warning('No overlapping ranges, returning NULL');
    return(NULL);
  }
  ov = as.data.frame(ov); 
  ov = split(ov, ov$queryHits);
  out = as.list(rep(NA, length(qGR)));
  qDF = as.data.frame(qGR);
  names(out) = paste0(qDF$seqnames, ':', qDF$start, '-', qDF$end);
  for (i in 1:length(ov)) {
    qind = as.numeric(names(ov)[i]);
    qr = as.data.frame(ranges(qGR[qind]));
    qname = paste0(as.character(seqnames(qGR[qind])), ':', qr$start, '-', qr$end);
    out[[qname]] = as.data.frame(sGR[ov[[i]]$subjectHits]);
  }
  return(out);
}

.getOverlapsStranded = function (qGR, sGR, ...) {
  ov = findOverlaps(qGR, sGR, maxgap=0L, ignore.strand=FALSE, ...);
  if (length(ov) == 0) {
    warning('No overlapping ranges, returning NULL');
    return(NULL);
  }
  ov = as.data.frame(ov); 
  ov = split(ov, ov$queryHits);
  out = as.list(rep(NA, length(qGR)));
  qDF = as.data.frame(qGR);
  names(out) = paste0(qDF$seqnames, ':', qDF$start, '-', qDF$end);
  for (i in 1:length(ov)) {
    qind = as.numeric(names(ov)[i]);
    qr = as.data.frame(ranges(qGR[qind]));
    qname = paste0(as.character(seqnames(qGR[qind])), ':', qr$start, '-', qr$end);
    out[[qname]] = as.data.frame(sGR[ov[[i]]$subjectHits]);
  }
  return(out);
}

# assumes that intervals are already known to overlap
.computeOverlapAmountUnstranded = function (s1, s2, e1, e2) {
  pos = sort(c(s1, s2, e1, e2));
  return(pos[3]-pos[2]+1);
}

.getOverlapsWithAmountsUnstranded = function (qGR, sGR, ...) {
  cat('Getting overlaps...\n');
  ovlist = .getOverlapsUnstranded(qGR, sGR, ...);
  loopinds = which(!sapply(ovlist, function(f) all(is.na(f))));
  for (i in loopinds) {
    thisdf = ovlist[[i]];
    thisdf = cbind(thisdf, ov.amount=rep(NA, nrow(thisdf)));
    qint = as.numeric(unlist(strsplit(unlist(strsplit(names(ovlist)[i], ':'))[2], '-')));
    for (row in 1:nrow(thisdf)) {
      thisdf$ov.amount[row] = .computeOverlapAmountUnstranded(qint[1], qint[2], 
                                                              thisdf$start[row], thisdf$end[row]);
    }
    ovlist[[i]] = thisdf;
  }
  return(ovlist);
}

# qrange should be one range in a GRanges object
.findNearestDownstreamFeatureUnstranded = function (qrange, sGR) {
  qdf = as.data.frame(qrange);
  qs = qdf$start; qe = qdf$end;
  sGR = sGR[ as.character(seqnames(sGR)) == as.character(qdf$seqnames) ];
  df = as.data.frame(cbind(as.data.frame(sGR), dist=distance(qrange, sGR)));
  df = df[order(df$dist), ];
  for (i in 1:nrow(df)) {
    if (df$strand[i] == '-') {
      if (qs > df$end[i]) {
        return(as.data.frame(cbind(df[i, ], skipped=(i-1))));
      }
    }
    if (df$strand[i] == '+') {
      if (qe < df$start[i]) {
        return(as.data.frame(cbind(df[i, ], skipped=(i-1))));
      }
    }
  }
  warning('No downstream neighbor for query, returning NA');
  return(NA);
}

.getNeighborsUnstranded = function (qGR, sGR, verbose=100) {
  cat('Getting nearest neighbors...\n');
  nn = .getNearestNeighborsUnstrandedNoOverlap(qGR, sGR);
  cat('second pass over queries that were downstream of both neighbors...\n');
  for (i in 1:length(nn)) {
    if (i %% verbose == 0) { cat(i,'...') }
    thisdf = nn[[i]];#print(thisdf)
    if (all(is.na(thisdf))) { next } 
    uprows = which(thisdf$qstream=='up');#print(uprows)
    if (length(uprows) > 0 && !is.na(thisdf$dist[uprows])) { next }
    qspl = unlist(strsplit(names(nn)[i], ':'));
    qint = as.numeric(unlist(strsplit(qspl[2], '-')));
    qrange = GRanges(seqnames=qspl[1], ranges=IRanges(start=qint[1], end=qint[2]), strand='*');
    nup = suppressWarnings(.findNearestDownstreamFeatureUnstranded(qrange, sGR));
    if (all(is.na(nup))) {
      warning(paste0('No downstream neighbor for ', names(nn)[i]));
      next;
    }
    tobind = c(nup[, c(6,1:3,5)], nup$dist, 'up');
    names(tobind) = names(thisdf);
    nn[[i]] = as.data.frame(rbind(thisdf, tobind));
    rownames(nn[[i]])[nrow(nn[[i]])] = 'nearest.down';
  }
  return(nn);
}

.getNAListElements = function (xlist) {
  return(sapply(xlist, function(f) all(is.na(f))));
}

# assumes mcols(GR) exists and has feature names
.getFeatureStartsFromGR = function (GR, mcolCol=1, stranded=FALSE) {
  starts = rep(NA, length(GR));
  strandvec = as.character(strand(GR));
  fwdrows = which(strandvec=='+');
  revrows = which(strandvec=='-');
  r = as.data.frame(ranges(GR));
  starts[fwdrows] = r$start[fwdrows];
  starts[revrows] = r$end[revrows];
  names(starts) = mcols(GR)[,mcolCol];
  strand = '*';
  if (stranded==TRUE) { strand = strand(GR) }
  gr = GRanges(seqnames=seqnames(GR), ranges=IRanges(start=starts, end=starts), strand=strand);
  mcols(gr) = mcols(GR);
  return(list(vec=starts, gr=gr));
}

.checkDmrExprFoldChangeAnticorrelated= function (dmrs, exprfcCol, dmrfcCol) {
  genefc = dmrs[, exprfcCol];
  dmrfc = dmrs[, dmrfcCol];
  expected = (genefc<0 & dmrfc>0) | (genefc>0 & dmrfc<0);
  return(expected);
}


.removeOverlappingRanges = function (gr1, gr2) {
  oldnames = names(gr1);
  # get overlaps
  fov = findOverlaps(gr1, gr2);
  # get actual overlapping ranges
  fovranges = ranges(fov, ranges(gr1), ranges(gr2));
  # split into lists to deal with queries that hit >1 subject
  fovspl = split(as.data.frame(fov), as.data.frame(fov)$queryHits);
  fovrangesspl = split(as.data.frame(fovranges), as.data.frame(fov)$queryHits);
  # loop through query hits
  contained = c();
  for (r in 1:length(fovspl)) {
    cat(unique(fovspl[[r]]$queryHits),'...')
    # get query and reduced subject range(s)
    qr = gr1[unique(fovspl[[r]]$queryHits)];
    sr = reduce(gr2[fovspl[[r]]$subjectHits]);
    # split into disjoined ranges
    qsdis = disjoin(c(qr, sr, ignore.mcols=T));
    # get this overlapping range
    rr = reduce(GRanges(seqnames=seqnames(qr), 
                        ranges=IRanges(start=fovrangesspl[[r]]$start, 
                                       end=fovrangesspl[[r]]$end)));
    qsdis = qsdis[-as.data.frame(findOverlaps(qsdis,rr))$queryHits];
    qsdis = qsdis[as.data.frame(findOverlaps(qsdis,qr))$queryHits];#print(as.data.frame(qsdis))
    if (length(qsdis) > 1) {
      qsdis = qsdis[order(-as.data.frame(qsdis)$width)][1];
      contained = c(contained, unique(fovspl[[r]]$queryHits))
    }
    ranges(gr1[unique(fovspl[[r]]$queryHits)]) <- ranges(qsdis);
  }
  x = as.data.frame(gr1);
  newnames = paste0(x$seqnames,':',x$start,'-',x$end);
  names(gr1) <- newnames;
  return(list(gr=gr1,oldnames=oldnames,ov=fov,contained=contained));
}



# dir shouldn't contain anything besides one subdir for each subject
.loadTophat2junctions = function (dir) {
  subdirs = list.files(dir);
  df = list();
  gr = list();
  for (f in 1:length(subdirs)) {
    cat(f,'...')
    df[[subdirs[f]]] = read.table(paste0(dir, subdirs[f], '/junctions.bed'), skip=1, sep='\t');
    gr[[subdirs[f]]] = GRanges(seqnames=df[[subdirs[f]]][,1], 
                               ranges=IRanges(start=df[[subdirs[f]]][,2], end=df[[subdirs[f]]][,3]), 
                               strand=df[[subdirs[f]]][,6]);
  }
  return(list(df=df,gr=gr))
}

.combineHandAndMachineAnnos = function (recipBL, handannos) {
  checkoverlap = handannos$gene %in% recipBL$gene;
  haInRecip = which(checkoverlap);
  haNoRecip = which(!checkoverlap);
  for (i in haInRecip) {
    reciprow = which(recipBL$gene == handannos$gene[i]);
    hsa1 = unlist(strsplit(recipBL$hsaHomologEntrez[reciprow], ','));
    hsa2 = unlist(strsplit(handannos$hsaHomologEntrez[i], ','));
    hsavec = unique(as.vector(na.omit(c(hsa1, hsa2))));
    recipBL$hsaHomologEntrez[reciprow] = paste0(hsavec, collapse=',');
  }
  newha = as.data.frame(matrix(nrow=length(haNoRecip), ncol=ncol(recipBL)));
  names(newha) = names(recipBL);
  newha$gene = handannos$gene[haNoRecip];
  newha$hsaHomologEntrez = handannos$hsaHomologEntrez[haNoRecip];
  newha$oBlastxHitDescription = handannos$description[haNoRecip];
  newrecipBL = as.data.frame(rbind(recipBL,newha));
  return(list(new=newrecipBL, haInRecip=haInRecip, haNoRecip=haNoRecip));
}

.getWindowsAroundGeneStarts = function (gr, up=5000, down=1000, sortOutput=FALSE) {
  if (sortOutput) {
    seqlevels(gr) <- sort(seqlevels(gr));
    gr = sort(gr);
  }
  df = as.data.frame(gr);
  df = df[, -which(names(df)=='width')];
  fwdrows = df$strand == '+';
  revrows = df$strand == '-';
  df$end[fwdrows] = df$start[fwdrows] + down;
  df$start[fwdrows] = df$start[fwdrows] - up;
  df$start[revrows] = df$end[revrows] - down;
  df$end[revrows] = df$end[revrows] + up;
  firstmcol = which(grepl('^mcols', names(df)))[1];
  names(df) = gsub('mcols.','',names(df),fixed=T);
  newmcols = df[, c(firstmcol:ncol(df))];
  newgr = GRanges(seqnames=df$seqnames, 
                  ranges=IRanges(start=df$start, end=df$end), 
                  strand=df$strand, 
                  mcols=newmcols);
  return(list(gr=newgr, df=as.data.frame(newgr)));
}

.fixInvalidPositions = function (gr, chrLengthVec) {
  df = as.data.frame(gr);
  df = as.data.frame(cbind(df, chrLength=chrLengthVec[match(df$seqnames, names(chrLengthVec))]));
  df$start[df$start < 0] = 1;
  df$end[df$end > df$chrLength] = df$chrLength[df$end > df$chrLength];
  firstmcol = which(grepl('^mcols', names(df)))[1];
  names(df) = gsub('mcols.','',names(df),fixed=T);
  newmcols = df[, c(firstmcol:ncol(df))];
  newgr = GRanges(seqnames=df$seqnames, 
                  ranges=IRanges(start=df$start, end=df$end), 
                  strand=df$strand, 
                  mcols=newmcols);
  return(list(gr=newgr, df=as.data.frame(newgr)));
}

.getFeatureFromNestedDmrHitList = function (dmrsL, featureName, hitType=c('ov','nn')) {
  hitType = hitType[1];
  if (!(hitType %in% c('ov','nn'))) { stop(paste0('hitType must be be \'ov\' or \'nn\'')) }
  fxn = paste0('f$',hitType,'$',featureName);
  fL = lapply(dmrsL, function(f) eval(parse(text=fxn)));
  return(fL);
}

.removeNaRowsFromDFsInList = function (listOfDataframes) {
  nrowDfs = lapply(listOfDataframes, nrow);
  inds = (1:length(listOfDataframes))[!sapply(nrowDfs, is.null)];
  for (i in inds) {
    removeRow = apply(listOfDataframes[[i]], 1, function(f) all(is.na(f)));
    listOfDataframes[[i]] = listOfDataframes[[i]][!removeRow, ];
  }
  return(listOfDataframes);
}

.subsetNestedNNHitListOnStreamAndDistance = function (dmrsGRhitsNN, stream='up', dist=20000) {
  x = lapply(lapply(dmrsGRhitsNN, 
                    function(ff) lapply(ff[!is.na(ff)], 
                                        function(f) f[which(f$qstream==stream & f$dist <= dist), ])), 
             function (fff) fff[sapply(fff,nrow) > 0]);
  return(x)
}

.searchForGenesAcrossNestedDmrHitList = function (dmrsGRhitsNN, genes, geneCol) {
  x = list()
  for (geneOI in genes) {
    y = lapply(lapply(dmrsGRhitsNN, 
                      function(ff) lapply(ff[!is.na(ff)], 
                                          function(f) f[which(f[,geneCol]==geneOI), ])), 
               function (fff) fff[sapply(fff,nrow) > 0]);
    x[[geneOI]] = y[sapply(y,length)>0];
    cat(paste0(geneOI, ': ', length(x[[geneOI]]), '\n'));
  }
  x = x[sapply(x, length) > 0];
  return(x);
}


.getDmrFeatureOverlapsAcrossGroups = function (dmrsGR, features, groups) {
  qname = deparse(substitute(dmrsGR));
  if (is.null(names(groups))) { 
    names(groups) = paste0('group', 1:length(groups));
  }
  for (feat in features) {
    cat(paste0('\n------------------------------------ ',feat,'...\n'));
    # nearest neighbors
    thisname = paste0(qname,'_',feat,'NN');
    assign(thisname, .getNeighborsUnstranded(dmrsGR, get(paste0(feat,'GR'))), pos=1);
    tmp = get(thisname, pos=1);
    for (g in 1:length(groups)) {
      assign(paste0(thisname,'.',names(groups)[g]), tmp[names(tmp) %in% groups[[g]]], pos=1);
    }
    # overlaps
    thisname = paste0(qname,'_',feat,'OV');
    assign(thisname, .getOverlapsWithAmountsUnstranded(dmrsGR, get(paste0(feat,'GR'))), pos=1);
    tmp = get(thisname, pos=1);
    for (g in 1:length(groups)) {
      assign(paste0(thisname,'.',names(groups)[g]), tmp[names(tmp) %in% groups[[g]]], pos=1);
    }
  }
}

.getDmrFeatureOverlapsAcrossGroupsNoNN = function (dmrsGR, features, groups) {
  qname = deparse(substitute(dmrsGR));
  if (is.null(names(groups))) { 
    names(groups) = paste0('group', 1:length(groups));
  }
  for (feat in features) {
    # overlaps
    thisname = paste0(qname,'_',feat,'OV');
    assign(thisname, .getOverlapsWithAmountsUnstranded(dmrsGR, get(paste0(feat,'GR'))), pos=1);
    tmp = get(thisname, pos=1);
    for (g in 1:length(groups)) {
      assign(paste0(thisname,'.',names(groups)[g]), tmp[names(tmp) %in% groups[[g]]], pos=1);
    }
  }
}

.getFeaturesFromDmrOverlapOrNeighborList = function (dmrsOV_NN, featureName, nestedFeatures=NULL) {
#   if (is.character(nestedFeatures)) {
#     fOV_NN = .getFeatureFromDmrOverlapHitList(dmrsOV_NN, featureName, nestedFeatures);
#   } else {
#     fOV_NN = dmrsOV_NN;
#   }
#  fOV_NN = dmrsOV_NN;
  fOV_NN = .removeNaRowsFromDFsInList(dmrsOV_NN);
  numi = lapply(fOV_NN, nrow);
  numi[which(sapply(numi, is.null))] = 0;
  numi = unlist(numi);
  nai = which(numi == 0);
  mi = which(numi > 1);
  if (length(nai) == 0) { 
    fvecinds = 1:length(fOV_NN);
  } else {
    fvecinds = -nai;
  }
  fvec = unlist(lapply(fOV_NN[fvecinds], function(f) f[, which(names(f)==featureName)]));
  return(list(fOV_NN=fOV_NN, fvec=fvec, numHits=numi, noHit=nai, multiHit=mi));
}

.subsetNearestNeighborHitListOnDistance = function (dmrsNN, up=NULL, down=NULL) {
  numi = lapply(dmrsNN, nrow);
  goodinds = which(!sapply(numi, is.null));
  for (i in goodinds) {
    keeprows = c();
    this = dmrsNN[[i]];
    if (is.numeric(up)) {
      keeprows = c(keeprows, which((this$qstream == 'up') & (this$dist <= up)));
    }
    if (is.numeric(down)) {
      keeprows = c(keeprows, which((this$qstream == 'down') & (this$dist <= down)));
    }
    dmrsNN[[i]] = this[keeprows, ];
  }
  return(dmrsNN);
}


# only works for genes
.combineOverlapAndNearestNeighborHitsForDmr = function (nnDF, ovDF, fnn='name', fov='mcols.geneSym') {
  nnr = nrow(nnDF); 
  ovr = nrow(ovDF);
  new = as.data.frame(matrix(nrow=1, ncol=12));
  names(new) = c('seqnames','start','end','strand','feature',
                 'dist','qstream',
                 'mcols.geneSym','mcols.geneID','mcols.biotype','mcols.GC','mcols.isoforms');
  if (!is.null(ovr)) {
    tobind = as.data.frame(cbind(seqnames=as.character(ovDF$seqnames),
                                 start=ovDF$start, end=ovDF$end, strand=as.character(ovDF$strand),
                                 feature=ovDF[, which(names(ovDF)==fov)],
                                 dist=NA, qstream=NA,
                                 mcols.geneSym=ovDF$mcols.geneSym, mcols.geneID=ovDF$mcols.geneID,
                                 mcols.biotype=ovDF$mcols.biotype, mcols.GC=ovDF$mcols.GC,
                                 mcols.isoforms=ovDF$mcols.isoforms));
    new = as.data.frame(rbind(new, tobind));
  } 
  if (!is.null(nnr)) {
    tobind = as.data.frame(cbind(seqnames=as.character(nnDF$seqnames),
                                 start=nnDF$start, end=nnDF$end, strand=as.character(nnDF$strand),
                                 feature=nnDF[, which(names(nnDF)==fnn)],
                                 dist=nnDF$dist, qstream=nnDF$qstream,
                                 mcols.geneSym=NA,mcols.geneID=NA,mcols.biotype=NA,mcols.GC=NA,mcols.isoforms=NA));
    new = as.data.frame(rbind(new, tobind));
  }

  return(new[-1,])
}
.combineOverlapAndNearestNeighborHits = function (dmrsNN, dmrsOV, fnn='name', fov='mcols.geneSym') {
  if (!all(names(dmrsNN)==names(dmrsOV))) { stop() }
  dmrnames = names(dmrsNN);
  new = list();
  for (i in 1:length(dmrnames)) {
    new[[dmrnames[i]]] = .combineOverlapAndNearestNeighborHitsForDmr(dmrsNN[[i]], dmrsOV[[i]], fnn, fov);
  }
  return(new);
}

# #features = c('gene', 'tss', 'br', 'mrna', 'mrnaSs', 'exon', 'cds', 'te', 'mi', 'z');
# #featureListNamesVec = paste0(features,'OV');
# .combineAllFeatureHitListsByDmr = function (featureListNamesVec) {
#   noexist = c();
#   for (f in featureListNamesVec) {if (!exists(f, mode='list')) { noexist = c(noexist, f) }}
#   if (length(noexist) > 0) {
#     fmiss = which(featureListNamesVec %in% noexist);
#     warning(paste0(paste0(featureListNamesVec[fmiss], collapse=', '),' not in the workspace'));
#     featureListNamesVec = featureListNamesVec[-fmiss];
#   }
#   hits = list();
#   dmrnames = names(get(featureListNamesVec[1]));
#   for (dmr in dmrnames) {
#     hits[[dmr]] = list(nn=list(),ov=list());
#     for (feat in featureListNamesVec) {
#       hits[[dmr]]
#     }
#   }
#   
# }

.limitVectorsInList = function (veclist, num=3, method=c('first','first.sorted','random'), smartUnlist=FALSE) {
  method = method[1];
  if (!(method %in% c('first','first.sorted','random'))) { stop() }
  tolimit = which(sapply(veclist, length) > num);
  if (method=='first') {
    veclist[tolimit] = lapply(veclist[tolimit], function(f) f[1:num]);
  } else if (method=='first.sorted') {
    veclist[tolimit] = lapply(veclist[tolimit], function(f) sort(f)[1:num]);
  } else {
    veclist[tolimit] = lapply(veclist[tolimit], function(f) sample(f, num));
  }
  outlist = list(new=veclist, limited=tolimit);
  if (smartUnlist) {
    outlist[['newvec']] = .smartUnlistVecListToVec(veclist, sep='__', preserveOrder=FALSE);
  }
  return(outlist);
}

# ~9-10x slower if preserveOrder=TRUE
.smartUnlistVecListToVec = function (veclist, sep='__', preserveOrder=FALSE) {
  if (preserveOrder) {
    newvec = c();
    for (v in 1:length(veclist)) {
      thisvec = veclist[[v]];
      if (length(thisvec) == 1) {
        newvec = c(newvec, thisvec);
        names(newvec)[length(newvec)] = names(veclist)[v];
      } else if (length(thisvec) > 1) {
        names(thisvec) = paste0(names(veclist)[v], sep, 1:length(thisvec));
        newvec = c(newvec, thisvec);
      } else {
        stop()
      }
    }
  } else {
    mi = which(sapply(veclist, length) > 1);
    if (length(mi) > 0) {
      newvec = unlist(veclist[-mi]);
      mlist = veclist[mi];
      for (mv in 1:length(mlist)) {
        thisvec = mlist[[mv]];
        names(thisvec) = paste0(names(mlist)[mv], sep, 1:length(thisvec));
        newvec = c(newvec, thisvec);
      }
    } else {
      newvec = unlist(veclist);
    }
  }
  return(newvec);
}

.buildIdsFromOV_NN = function (fOV_NN, featureName, maplist, limit, limitMethod, smartUnlist) {
  newfOV_NN = .getFeaturesFromDmrOverlapOrNeighborList(dmrsOV_NN=fOV_NN, 
                                                       featureName=featureName,
                                                       nestedFeatures=NULL);
  ids0 = .mapAbGenesHsEntrez(ids=newfOV_NN$fvec, maplist=maplist, inputType='ab');
  ids = .limitVectorsInList(veclist=ids0, num=limit, method=limitMethod, smartUnlist=smartUnlist);
  return(list(processedInput=newfOV_NN, ids=ids));
}

# 
# .buildGeneListsForGoFromOvOrNn = function (fOV_NN, extractFeatures, groups=NULL, maplist=NULL, limit=NULL, limitMethod=NULL) {
#   l = .getFeaturesFromDmrOverlapOrNeighborList(fOV_NN, extractFeatures, NULL);
#   ids = list(all=l$fvec);
#   if (is.list(groups)) {
#     if (is.null(names(groups))) { names(groups) = paste0('group',1:length(groups)) }
#     for (g in 1:length(groups)) {
#       tmp = l$fOV_NN[names(l$fOV_NN) %in% groups[[g]]];
#       ids[[names(groups)[g]]] = .getFeaturesFromDmrOverlapOrNeighborList(tmp, extractFeatures, NULL)$fvec;
#     }
#   }
#   if (is.list(maplist)) {
#     if (is.null(limit)) { limit=2 }
#     if (is.null(limitMethod)) { limitMethod='first.sorted' }
#     idsmap = lapply(ids, function(f) .mapAbGenesHsEntrez(unique(f), maplist, 'ab'));
#     idsmap = lapply(idsmap, function(f) .limitVectorsInList(f, limit, limitMethod));
#     ids = list(ids=ids, idsmap=idsmap);
#   }
#   return(ids);
# }
# 

.getSymmetricWindowsAroundFeatures = function (gr, bpToAdd, chrLengthVec, sortOutput=FALSE) {
  if (sortOutput) { gr = .sortGRanges(gr) }
  r = as.data.frame(ranges(gr));
  newgr = GRanges(seqnames=seqnames(gr),
                  ranges=IRanges(start=(r$start - bpToAdd),
                                 end=(r$end = r$end + bpToAdd)),
                  strand=strand(gr),
                  mcols=mcols(gr));
  newgr = .fixInvalidPositions(newgr, chrLengthVec)$gr;
  return(newgr)
}


.buildTfHitSummary = function (tfhitsDmrsMulti) {
  dflist = lapply(tfhitsDmrsMulti, function(f) f$hits$df);
  outnames = c('ID','TF','class','length','hits','mean.absScore','mean.relScore');
  out = as.data.frame(matrix(nrow=length(dflist),ncol=length(outnames)));
  names(out) = outnames;
  for (i in 1:length(dflist)) {
    this = dflist[[i]];
    if (is.null(this)) { next }
    out$ID[i] = this$ID[1];
    out$TF[i] = this$TF[1];
    out$class[i] = this$class[1];
    out$length[i] = ncol(attributes(tfhitsDmrsMulti[[i]]$tfm$pfm)$profileMatrix);
    out$hits[i] = nrow(this);
    out$mean.absScore[i] = mean(this$absScore);
    out$mean.relScore[i] = mean(this$relScore);
  }
  return(out[!apply(out,1,function(f) all(is.na(f))),]);
}

.groupMultiTfbsHitsByDmrs = function (tfhitsDmrsMulti, orderHits=TRUE, verbose=100) {
  numdmrs = unique(sapply(tfhitsDmrsMulti, function(f) length(f$hits$raw)));
  if (length(numdmrs) != 1) { stop() }
  dmrlist = vector(mode='list',length=numdmrs);
  names(dmrlist) = names(tfhitsDmrsMulti[[1]]$hits$raw);
  if (orderHits) {
    for (d in 1:numdmrs) {
      if (d %% verbose == 0) { cat(d, '...')}
      dmrlist[[d]] = do.call('rbind', lapply(tfhitsDmrsMulti, function(f) f$hits$raw[[d]]));
      if (nrow(dmrlist[[d]]) > 0) { dmrlist[[d]] = dmrlist[[d]][order(-dmrlist[[d]]$absScore), ] }
    }
  } else {
    for (d in 1:numdmrs) {
      if (d %% verbose == 0) { cat(d, '...')}
      dmrlist[[d]] = do.call('rbind', lapply(tfhitsDmrsMulti, function(f) f$hits$raw[[d]]));
    }
  }
  return(dmrlist);
}

.getLociInDmr = function (chr, start, end, BSseq, seq=NULL, ...) {
  if (class(BSseq) == 'BSseq') {
    df = as.data.frame(attributes(BSseq[seqnames(BSseq)==chr])$rowData);
  } else if (class(BSseq) == 'BSseqTstat') {
    df = as.data.frame(attributes(BSseq)$gr);
  } else {
    stop();
  }
  abspos = df$start[df$start >= start & df$end <= end];
  dmrpos = abspos - start + 1;
  nucs = c();
  if (class(seq) == 'DNAString') {
    for (i in dmrpos) { 
      nucs = c(nucs, as.character(subseq(seq, i, i)));
    }
  }
  return(list(abspos=abspos,dmrpos=dmrpos,nucs=nucs));
}

.countLociInTfHit = function (tfHitsForDmr, chr, start, end, BSseq, ...) {
  dmrloci = .getLociInDmr(chr=chr, start=start, end=end, BSseq=BSseq, ...);
  dmrlociGR = GRanges(seqnames=chr, ranges=IRanges(start=dmrloci$dmrpos, end=dmrloci$dmrpos));
  tfhitGR = GRanges(seqnames=chr, ranges=IRanges(start=tfHitsForDmr$start, end=tfHitsForDmr$end));
  ov = .getOverlapsUnstranded(tfhitGR, dmrlociGR);
  ov = ov[!sapply(ov, function(f) all(is.na(f)))];
  numHits = sapply(ov, nrow);
  return(list(dmrloci=dmrloci, ov=list(hits=ov,numHits=numHits)));
}

.getTfHitOverlaps = function (tfHitsForDmr, ...) {
  gr = GRanges(seqnames=unique(sapply(strsplit(tfHitsForDmr$seqnames,':'), function(f) f[1])),
               ranges=IRanges(start=tfHitsForDmr$start, 
                              end=tfHitsForDmr$end),
               mcols=data.frame(absScore=tfHitsForDmr$absScore,
                                relScore=tfHitsForDmr$relScore,
                                ID=tfHitsForDmr$ID, 
                                TF=tfHitsForDmr$TF));
  ov = .getOverlapsWithAmountsUnstranded(gr, gr, ...);
  names(ov) = paste0(names(ov),'_',mcols(gr)$mcols.TF);
  return(list(tfHitGR=gr, tfHitOV=ov));
}

.buildTfHitOverlapTable = function (getTfHitOverlapsOutput, motifLength=NULL, thresh=NULL) {
  l0 = getTfHitOverlapsOutput;
  #names(l0[[2]]) = paste0(names(l0[[2]]),'_',mcols(l0[[1]])$mcols.TF);
  l = l0[[2]];
  n = length(l);
  df = as.data.frame(matrix(rep(0,n^2), ncol=n, nrow=n,
                            dimnames=list(names(l),names(l))));
  naids = c();
  for (h in 1:length(l)) {
    this = l[[h]];#print(h)
    if (!all(is.na(this))) {
      #print(this)
      ids = (paste0(this$seqnames,':',this$start,'-',this$end,'_',this$mcols.TF)); #print(ids)
      if (any(duplicated(ids))) {#print(ids)
        this = this[!duplicated(ids), ];#print(this)
        ids = ids[!duplicated(ids)];#print(ids)
      }
      #df[h, match(ids, names(df))] = 1;
      # print(match(ids, names(df)))
      #   print(this)
      if (is.numeric(motifLength)) {
        keep = which(this$width >= motifLength);
        this = this[keep, ];
        ids = ids[keep];
      }
      df[h, match(ids, names(df))] = this$ov.amount/this$width;
      df[match(ids, names(df)), h] = this$ov.amount/this$width;
    } else {
      naids = c(naids,h);
    }
  }
  df = df[-h,-h];
  if (is.numeric(thresh)) {
    df[which(df < thresh, arr.ind=T)] = 0;
  }
  return(df);
}

.buildTfIdNameLookup = function (tfhitsDmrsMulti) {
  return(sapply(tfhitsDmrsMulti, function(f) attributes(f$tfm[[1]])$name));
}

.countLociInTfHitAcrossDmrs = function (groupMultiTfbsHitsByDmrsOutput, BSseq, verbose=100) {
  dl = groupMultiTfbsHitsByDmrsOutput;
  l = vector(mode='list',length=length(dl));
  names(l) = names(dl);
  loopinds = which(sapply(dl, nrow) > 0);
  for (i in loopinds) {
    if (i %% verbose == 0) { cat(i,'...') }
    thisname = strsplit(names(dl)[i], ':')[[1]];
    thesecoords = as.numeric(strsplit(thisname[2],'-')[[1]]);
    l[[names(dl)[i]]] = .countLociInTfHit(tfHitsForDmr=dl[[i]],
                                          chr=thisname[1], 
                                          start=thesecoords[1], 
                                          end=thesecoords[2], 
                                          BSseq=BSseq);
  }
  return(l);
}

.analyzeTfHitsByDmr = function (tfhitsDmrsMulti, BSseq, verbose=100) {
  stats = .buildTfHitSummary(tfhitsDmrsMulti);
  cat('-------------------- Grouping hits by dmr...\n');
  gByDmr = .groupMultiTfbsHitsByDmrs(tfhitsDmrsMulti, verbose=verbose);
  nHitsByDmr = sapply(gByDmr, nrow);
  cat('\n\n-------------------- Counting methylation loci in each tfbs hit...\n');
  cpgByHit = .countLociInTfHitAcrossDmrs(gByDmr, BSseq, verbose=verbose);
  maxcpg = rep(NA, length(cpgByHit));
  names(maxcpg) = names(cpgByHit);
  loopinds = which(!(sapply(cpgByHit, function(f) is.null(names(f)))));
  for (i in loopinds) {
    if (class(cpgByHit[[i]]$ov$numHits) == 'list') {next}
    maxcpg[i] = max(cpgByHit[[i]]$ov$numHits);
  }
  gByDmr = .addNumLociToTfHitDfForDmr(gByDmr, cpgByHit);
  cat('\n\n-------------------- Collecting tfbs hits that overlap cpgs in dmrs...\n');
  tfcpg = vector(mode='list',length=length(gByDmr));
  names(tfcpg) = names(gByDmr);
  loopinds = which(sapply(gByDmr, function(f) !all(is.null(f$cpgs[[1]]))));
  for (i in loopinds) {
    if (i %% verbose == 0) { cat(i,'...') }
    tfcpg[[i]] = gByDmr[[i]][which(gByDmr[[i]]$cpg > 0), ];
  }
  return(list(stats=stats, byDmr=gByDmr, nHits=nHitsByDmr, cpgByHit=cpgByHit, maxCpgByHit=maxcpg, hitsWithCpGs=tfcpg));
}







.addNumLociToTfHitDfForDmr = function (tfhits, numloci) {
  if (any(names(tfhits) != names(numloci))) { stop() }
  loopinds = which(sapply(tfhits,nrow) != 0);
  for (i in loopinds) {
    this = tfhits[[i]];
    this = as.data.frame(cbind(this, cpgs=rep(0,nrow(this))));
    thischr = strsplit(this$seqnames, ':')[[1]][1];
    hitnames = paste0(thischr, ':', this$start, '-', this$end);#print(hitnames)
    thisnumloci = numloci[[i]]$ov$numHits;#print(thisnumloci)
    this$cpgs = thisnumloci[match(hitnames, names(thisnumloci))];#print(this$cpgs)
    this$cpgs[is.null(this$cpgs)] = 0; #careful of NULLs
    tfhits[[i]] = this;
  }
  return(tfhits);
}

.getTstatsForPositionsOnChr = function (posvec, BSseqTstatChr) {
  tstatschr = .getStatsFromBSseqTstat(BSseqTstatChr);
  return(tstatschr[match(posvec, tstatschr$start), ]);
}

.getTstatsForRange = function (chr, start, end, BSseqTstatChr, ...) {
  tt = .getLociInDmr(chr=chr, start=start, end=end, BSseq=BSseqTstatChr, ...);
  return(.getTstatsForPositionsOnChr(tt$abspos, BSseqTstatChr));
}

.getTstatsForPositionsOnChrAcrossSmoothings = function (chr, posvec, BSseqTstatList, type=c('pos','range')) {
  lpos = sapply(BSseqTstatList, function(f) which(names(f) == chr));
  if (is.list(lpos)) { stop() }
  tl = list();
  type = type[1];
  if (type == 'pos') {
    for (i in 1:length(lpos)) {
      thist = BSseqTstatList[[i]][[lpos[i]]];
      tl[[names(BSseqTstatList)[i]]] = .getTstatsForPositionsOnChr(posvec, thist);
    }
  } else if (type == 'range') {
    if (length(posvec) > 2) { warning('using extreme values from posvec as range') }
    r1 = min(posvec); 
    r2 = max(posvec);
    for (i in 1:length(lpos)) {
      thist = BSseqTstatList[[i]][[lpos[i]]];
      # print(chr);print(r1);print(r2);print(thist)
      tl[[names(BSseqTstatList)[i]]] = .getTstatsForRange(chr, r1, r2, thist);
    }
  } else {
    stop('Invalid type argument')
  }
  avgt = vector(length=2, mode='numeric');
  avgt[1] = mean(sapply(tl, function(f) mean(f$tstat)));
  avgt[2] = mean(sapply(tl, function(f) mean(f$tstat.corrected)));
  names(avgt) = c('tstat','tstat.corrected');
  return(list(tstats=tl, means=avgt));
}







.verboseScatterplotAllColumnPairs = function (mat, mfrow=c(3,ncol(mat)), col=NULL, labels=NULL, ...) {
  if (is.null(col)) {
    col = rownames(mat);
  } else {
    col = col;
  }
  par(mfrow=mfrow);
  for (i in 1:ncol(mat)) {
    if (i==ncol(mat)) {break};
    for (j in (i+1):ncol(mat)) {
      WGCNA::verboseScatterplot(mat[,i], mat[,j], 
                         xlab=colnames(mat)[i], ylab=colnames(mat)[j],
                         frame.plot=F, col=col,
                         abline=T, ...
      );
      if (!is.null(labels)) {
        text(mat[,i], mat[,j], labels=labels)
      }
    }
  }
}

#.verboseScatterplotAllColumnPairsLabels = function (mat, mfrow=c(3,ncol(mat)), col=NULL, ...) {
#   if (is.null(col)) {
#     col = rownames(mat);
#   } else {
#     col = col;
#   }
#   par(mfrow=mfrow);
#   for (i in 1:ncol(mat)) {
#     if (i==ncol(mat)) {break};
#     for (j in (i+1):ncol(mat)) {
#       verboseScatterplot(mat[,i], mat[,j], 
#                          xlab=colnames(mat)[i], ylab=colnames(mat)[j],
#                          frame.plot=F, col=col,
#                          abline=T, ...
#       );
#     }
#   }
# }




.verboseScatterplotVecAgainstColumns = function (mat, vec, labels=NULL, xlab=NULL, ...) {
  if (any(names(vec) != rownames(mat))) {
    stop('names(vec) must be identical to rownames(mat)');
  }
  par(mfrow=c(2,ceiling(ncol(mat)/2)));
  for (i in 1:ncol(mat)) {
    if (is.null(xlab)) {
      xlab = deparse(substitute(vec))
    } else {
      xlab = xlab
    }
    WGCNA::verboseScatterplot(vec, mat[,i], 
                       xlab=xlab, ylab=colnames(mat)[i],
                       frame.plot=F, abline=T, ...
    );
    if (!is.null(labels)) {
      text(vec, mat[,i], labels=labels)
    }
  }
}




# intervals must be formatted as 'scaffold_0:1-100'
.write_bedfile_for_intervals = function (intervals, filename, returnDataFrame=T) {
  intsplit = strsplit(intervals, ':');
  starts_ends = sapply(intsplit, function(f) strsplit(f[2],'-'));
  df = as.data.frame(matrix(nrow=length(intervals), ncol=3));
  df[, 1] = sapply(intsplit, function(f) f[1]);
  df[, 2] = sapply(starts_ends, function(f) as.numeric(f[1]));
  df[, 3] = sapply(starts_ends, function(f) as.numeric(f[2]));
  if (!is.null(names(intervals))) {
    df = as.data.frame(cbind(df, names(intervals)));
  }
  write.table(df, file=filename, quote=F, col.names=F, row.names=F, sep='\t');
  if (returnDataFrame) {
    return(df);
  }
}


.buildGR_for_intervals = function (intervals) {
  intsplit = strsplit(intervals, ':');
  starts_ends = sapply(intsplit, function(f) strsplit(f[2],'-'));
  gr = GRanges(seqnames=sapply(intsplit, function(f) f[1]),
               ranges=IRanges(start=sapply(starts_ends, function(f) as.numeric(f[1])), 
                              end=sapply(starts_ends, function(f) as.numeric(f[2]))));
  return(gr)
}



.draw_pairwise_venn_from_groups = function (gp1, gp2, groupNames=NULL, ...) {
  gp1_2 = intersect(gp1, gp2);
  if (is.null(groupNames)) {
    groupNames = c(paste0(deparse(substitute(gp1)),':',length(gp1)), 
                   paste0(deparse(substitute(gp2)),':',length(gp2)));
  } else {
    for (g in 1:2) {
      groupNames[g] = paste0(groupNames[g],':',length(get(paste0('gp',g))));
    }
  }
  draw.pairwise.venn(length(gp1), length(gp2), length(gp1_2), 
                     category=groupNames,
                     fontfamily=rep('Helvetica',3),
                     cat.fontfamily=rep('Helvetica',2),
                     ...);
  uniques = list(setdiff(gp1, gp2), setdiff(gp2, gp1));
  groupNames = sapply(strsplit(groupNames,':'), function(f) f[1]);
  names(uniques) = groupNames;
  return(list(overlaps=gp1_2, uniques=uniques, groupNames=groupNames));
}

.draw_quad_venn_from_groups = function (gp1, gp2, gp3, gp4, groupNames=NULL, ...) {
  gp1_2 = intersect(gp1, gp2);
  gp1_3 = intersect(gp1, gp3);
  gp1_4 = intersect(gp1, gp4);
  gp2_3 = intersect(gp2, gp3);
  gp2_4 = intersect(gp2, gp4);
  gp3_4 = intersect(gp3, gp4);
  gp1_2_3 = intersect(gp1_2, gp3);
  gp1_2_4 = intersect(gp1_2, gp4);
  gp1_3_4 = intersect(gp1_3, gp4);
  gp2_3_4 = intersect(gp2_3, gp4);
  gp1_2_3_4 = intersect(gp1_2_3, gp4);
  if (is.null(groupNames)) {
    groupNames = c(paste0(deparse(substitute(gp1)),':',length(gp1)), 
                   paste0(deparse(substitute(gp2)),':',length(gp2)), 
                   paste0(deparse(substitute(gp3)),':',length(gp3)), 
                   paste0(deparse(substitute(gp4)),':',length(gp4)));
   # groupNames = gsub('\\)', '', gsub('rownames\\(', '', groupNames));
  } else {
    for (g in 1:4) {
      groupNames[g] = paste0(groupNames[g],':',length(get(paste0('gp',g))));
    }
  }
  draw.quad.venn(length(gp1), length(gp2), length(gp3), length(gp4),
                 length(gp1_2), length(gp1_3), length(gp1_4),
                 length(gp2_3), length(gp2_4),
                 length(gp3_4),
                 length(gp1_2_3), length(gp1_2_4), length(gp1_3_4),
                 length(gp2_3_4),
                 length(gp1_2_3_4),
                 category=groupNames,
                 fontfamily=rep('Helvetica',15),
                 cat.fontfamily=rep('Helvetica',4),
                 ...
  );
  overlaps = list(gp1_2, gp1_3, gp1_4, gp2_3, gp2_4, gp3_4, 
                  gp1_2_3, gp1_2_4, gp1_3_4, gp2_3_4, gp1_2_3_4);
  uniques = list(setdiff(gp1, c(gp2,gp3,gp4)), 
                 setdiff(gp2, c(gp1,gp3,gp4)),
                 setdiff(gp3, c(gp1,gp2,gp4)),
                 setdiff(gp4, c(gp1,gp2,gp3)));
  groupNames = sapply(strsplit(groupNames,':'), function(f) f[1]);
  names(uniques) = groupNames;
  names(overlaps) = c(paste0(groupNames[1:2],collapse=':'), 
                      paste0(groupNames[c(1,3)],collapse=':'), 
                      paste0(groupNames[c(1,4)],collapse=':'), 
                      paste0(groupNames[c(2,3)],collapse=':'),
                      paste0(groupNames[c(2,4)],collapse=':'),
                      paste0(groupNames[c(3,4)],collapse=':'),
                      paste0(groupNames[1:3],collapse=':'),
                      paste0(groupNames[c(1,2,4)],collapse=':'),
                      paste0(groupNames[c(1,3,4)],collapse=':'),
                      paste0(groupNames[2:4],collapse=':'),
                      paste0(groupNames[1:4],collapse=':')
  );
  return(list(overlaps=overlaps, uniques=uniques, groupNames=groupNames));
}





.draw_triple_venn_from_groups = function (gp1, gp2, gp3, groupNames=NULL, ...) {
  gp1_2 = intersect(gp1, gp2);
  gp1_3 = intersect(gp1, gp3);
  gp2_3 = intersect(gp2, gp3);
  gp1_2_3 = intersect(gp1_2, gp3);
  if (is.null(groupNames)) {
    groupNames = c(paste0(deparse(substitute(gp1)),':',length(gp1)), 
                   paste0(deparse(substitute(gp2)),':',length(gp2)), 
                   paste0(deparse(substitute(gp3)),':',length(gp3)));
   # groupNames = gsub('\\)', '', gsub('rownames\\(', '', groupNames));
  } else {
    for (g in 1:3) {
      groupNames[g] = paste0(groupNames[g],':',length(get(paste0('gp',g))));
    }
  }
  draw.triple.venn(length(gp1), length(gp2), length(gp3), 
                 length(gp1_2), length(gp1_3), 
                 length(gp2_3),
                 length(gp1_2_3),
                 category=groupNames,
                 fontfamily=rep('Helvetica',7),
                 cat.fontfamily=rep('Helvetica',3),
                 ...
  );
  overlaps = list(gp1_2, gp1_3, gp2_3, 
                  gp1_2_3);
  uniques = list(setdiff(gp1, c(gp2,gp3)), 
                 setdiff(gp2, c(gp1,gp3)),
                 setdiff(gp3, c(gp1,gp2)));
  groupNames = sapply(strsplit(groupNames,':'), function(f) f[1]);
  names(uniques) = groupNames;
  names(overlaps) = c(paste0(groupNames[1:2],collapse=':'), 
                      paste0(groupNames[c(1,3)],collapse=':'), 
                      paste0(groupNames[c(2,3)],collapse=':'),
                      paste0(groupNames[1:3],collapse=':')
  );
  return(list(overlaps=overlaps, uniques=uniques, groupNames=groupNames));
}



.screen_enrichGOAllOntologies = function (idlist, bgvec, OrgDb=org.Hs.eg.db, fdrmethvec=c('BY','BH'), fdrthvec=c(.01, .05)) {
  golist = list();
  for (fdrmeth in fdrmethvec) {
    cat(paste0('------------------------- ', fdrmeth,  ' -------------------------'));
    for (fdrth in fdrthvec) {
      cat(paste0('------------------------- ', fdrth,  ' -------------------------'));
      for (i in 1:length(idlist)) {
        cat(paste0('------------------------- ', names(idlist)[i],  ' -------------------------'));
        thisres = .enrichGOAllOntologies(idlist[[i]], OrgDb=OrgDb, refGenes=bgvec, pAdjustMethod=fdrmeth, pvalueCutoff=fdrth);
        thisname = paste0(fdrmeth,'.', gsub('0.','',fdrth,fixed=T));
        golist[[names(idlist)[i]]][[thisname]] = thisres;
      }
    }
  }
  return(golist);
}



.screen_enrichKEGG = function (idlist, bgvec, fdrmethvec=c('BY','BH','fdr'), fdrthvec=c(.01, .05, .1), ...) {
  kegglist = list();
  for (fdrmeth in fdrmethvec) {
    cat(paste0('\n------------------------- ', fdrmeth,  ' -------------------------'));
    for (fdrth in fdrthvec) {
      cat(paste0('\n  ', fdrth,  '...\n'));
      for (i in 1:length(idlist)) {
        cat(paste0('     ', names(idlist)[i], '...'));
        thisres = as.data.frame(summary(enrichKEGG(gene=idlist[[i]], universe=bgvec, pAdjustMethod=fdrmeth, pvalueCutoff=fdrth, ...)));
        thisname = paste0(fdrmeth,'.', gsub('0.','',fdrth,fixed=T));
        kegglist[[names(idlist)[i]]][[thisname]] = thisres;
      }
    }
  }
  return(kegglist);
}


.addDmrInfoToExpression = function (expr, dmrsOVbyGene) {
  alldmrgenes = unique(unlist(lapply(dmrsOVbyGene, names)));
  dmrcols = c('dmr.num', 'dmrfc', 'dmrfc.abs', 
              'n', 'width.dmr', 'invdensity', 
              'areaStat', 'areaStat.abs', 'meanDiff', 'meanDiff.abs');
  smat = matrix(ncol=length(dmrcols),nrow=nrow(expr),
                dimnames=list(rownames(expr),dmrcols));
  expr = as.data.frame(cbind(expr, smat));
  for (g in alldmrgenes) {
    row = match(g, rownames(expr));
    gdmr = lapply(dmrsOVbyGene, 
                  function(f) f[[match(g, names(f))]]);
    gdmr = gdmr[!sapply(gdmr, is.null)];
    gdmr = subset(do.call('rbind', gdmr), 
                  !duplicated(paste0(chr,start,end)));
    #
    expr$dmr.num[row] = nrow(gdmr);
    expr$dmrfc[row] = mean(log2(gdmr$group2.mean/gdmr$group1.mean));
    expr$dmrfc.abs[row] = mean(abs(log2(gdmr$group2.mean/gdmr$group1.mean)));
    expr$n[row] = mean(gdmr$n);
    expr$width.dmr[row] = mean(gdmr$width);
    expr$invdensity[row] = mean(gdmr$invdensity);
    expr$areaStat[row] = mean(gdmr$areaStat);
    expr$areaStat.abs[row] = mean(abs(gdmr$areaStat));
    expr$meanDiff[row] = mean(gdmr$meanDiff);
    expr$meanDiff.abs[row] = mean(abs(gdmr$meanDiff));
  }
  expr$dmr.num[is.na(expr$dmr.num)] = 0;
  fmat = matrix(nrow=nrow(expr), ncol=length(dmrsOVbyGene),
                dimnames=list(rownames(expr), names(dmrsOVbyGene)));
  for (f in 1:length(dmrsOVbyGene)) {
    fmat[, f] = rownames(expr) %in% names(dmrsOVbyGene[[f]]);
  }
  fmat = as.data.frame(cbind(fmat, 
                             fg=apply(fmat, 1, 
                                      function(f) paste0(colnames(fmat)[f],
                                                         collapse=':'))));
  fmat$fg[fmat$fg==''] = NA;
  expr = as.data.frame(cbind(expr, fmat));
  return(expr);
}


.filter_nullints_fromDMRs = function (dmrs, allnullints, overlapLimit=0.5, verbose=100) {
  nov = list();
  if (class(dmrs) != 'GRanges') {
    mcolCols = which(names(dmrs)=='n'):ncol(dmrs);
    dmrsGR = GRanges(seqnames=dmrs$chr, 
                     ranges=IRanges(start=dmrs$start, end=dmrs$end), 
                     mcols=dmrs[,mcolCols]);
    names(dmrsGR) = paste0(dmrs$chr, ':', dmrs$start, '-', dmrs$end);
  }
  if (class(allnullints) != 'GRanges') {
    allnullintsGR = .buildGR_for_intervals(allnullints);
  }
  for (d in 1:nrow(dmrs)) {
    if (d %% verbose == 0) { cat(d,'...') }
    thisw = dmrs$width[d];
    nov[[dmrsints[d]]] = as.data.frame(findOverlaps(dmrsGR[d], 
                                                    allnullintsGR, 
                                                    minoverlap=(thisw * overlapLimit)))
  }
  return(dmrs[!(names(dmrsGR) %in% names(nov)[sapply(nov, nrow)!=0]), ])
}








.getWindowsAroundGeneEnds = function (gr, down=5000, up=1000, sortOutput=FALSE) {
  ofwd = which(strand(gr)=='+');
  orev = which(strand(gr)=='-');
  strand(gr)[ofwd] <- '-';
  strand(gr)[orev] <- '+';
  newgr = .getWindowsAroundGeneStarts(gr=gr, up=down, down=up, sortOutput=sortOutput)$gr;
  strand(newgr)[ofwd] <- '+';
  strand(newgr)[orev] <- '-';
  return(list(gr=newgr, df=as.data.frame(newgr)));
}



.combineOvAndNnByDmr = function (features) {
  ov = paste0(features,'OV');
  nn = paste0(features,'NN');
  features = intersect(c(ov,nn), ls(pos=1));
  dmrints = names(get(features[1]));
  cat('Checking that all objects contain exact same dmrs...\n    ');
  for (feat in features) { 
    cat(feat,'...'); 
    if (any(names(get(feat)) != dmrints)) { stop() } 
  }
  dmrsGRhits = list();
  for (dmrint in dmrsints) {
    dmrsGRhits[[dmrint]] = list(nn=list(),ov=list());
    for (feat in grep('OV$',features,value=T)) {
      f = gsub('OV','',feat);
      dmrsGRhits[[dmrint]]$ov[[f]] = get(feat)[[dmrint]];
    }
    for (feat in grep('NN$',features,value=T)) {
      f = gsub('NN','',feat);
      dmrsGRhits[[dmrint]]$nn[[f]] = get(feat)[[dmrint]];
    }
  }
  return(dmrsGRhits);
}

.combineOvByDmr = function (features) {
  ov = paste0(features,'OV');
  features = intersect(c(ov), ls(pos=1));
  dmrints = names(get(features[1]));#print(head(dmrints))
  cat('Checking that all objects contain exact same dmrs...\n    ');
  for (feat in features) { 
    cat(feat,'...'); 
    if (any(names(get(feat)) != dmrints)) { stop() } 
  }
  dmrsGRhits = list();
  for (i in 1:length(dmrsints)) {#print('------------------------------------------------------')
    dmrsGRhits[[dmrints[i]]] = list();#print(dmrsGRhits[dmrint])
    for (feat in grep('OV$',features,value=T)) {
      f = gsub('OV','',feat);#print(f)
      
      tmp = get(feat);
      dmrsGRhits[[dmrints[i]]][[f]] = tmp[[i]];
      
    }
  }
  return(dmrsGRhits);
}

.make_numOV = function (dmrsGRhitsOV, features, prefix='dmrsGR_') {
  numOV = as.data.frame(matrix(nrow=nrow(dmrs), ncol=length(features)));
  dimnames(numOV) = list(dmrsints, features);
  for (feat in features) {
    featOV = lapply(dmrsGRhitsOV, function(f) f[[which(names(f) == paste0(prefix,feat))]] );
    featOVnums = sapply(featOV, nrow);
    featOVnums[which(sapply(featOVnums, is.null))] = 0;
    featOVnums = unlist(featOVnums);
    numOV[, match(feat, names(numOV))] = featOVnums;
  }
  return(numOV)
}


.generateIntervalsFromScaffold = function (scaffold_name, scaffold_length, interval_widths) {
  maxpos = scaffold_length - max(interval_widths);
  starts = sample(1:maxpos, length(interval_widths));
  ends = starts + interval_widths;
  gr = GRanges(seqnames=scaffold_name, ranges=IRanges(start=starts,end=ends));
  return(gr);
}

.process_dmrsGR_featureOV = function (dmrsGR_featureOV, dmrs, nameCol) {
  featureOVbyDmr = dmrsGR_featureOV[!is.na(dmrsGR_featureOV)];
  featureOVbyDmrSyms = lapply(featureOVbyDmr, function(f) unique(f[,nameCol])); 
  dmrOVbyFeature = do.call('rbind', featureOVbyDmr);
  dmrOVbyFeature = split(dmrOVbyFeature, dmrOVbyFeature[,nameCol]);
  # replace data frames with dmr info instead of gene info
  for (g in 1:length(dmrOVbyFeature)) {
    xdmrnames = strsplit(rownames(dmrOVbyFeature[[g]]), '.', fixed=T);
    xdmrnames = sapply(xdmrnames, function(f) f[1]);
    dmrOVbyFeature[[g]] = dmrs[rownames(dmrs) %in% xdmrnames, ];
  }
  return(list(byDmr=featureOVbyDmr, byDmrSyms=featureOVbyDmrSyms, byFeature=dmrOVbyFeature));
}



.writeBedsForTfHits = function (tfhits, stranded, chrCol, startCol, endCol, scoreCol, nameCol, filenameBase) {
  bedGraphCols = c(chrCol, startCol, endCol, scoreCol);
  bedCols = c(chrCol, startCol, endCol, nameCol, scoreCol);
  if (stranded) {
    write.table(tfhits[tfhits$strand=='+', bedGraphCols], 
                paste0(filenameBase,'_fwd.bedGraph'), quote=F, row.names=F, col.names=F, sep='\t');
    write.table(tfhits[tfhits$strand=='-', bedGraphCols], 
                paste0(filenameBase,'_rev.bedGraph'), quote=F, row.names=F, col.names=F, sep='\t');
  } else {
    write.table(tfhits[, bedGraphCols], 
                paste0(filenameBase,'.bedGraph'), quote=F, row.names=F, col.names=F, sep='\t');
  }
  write.table(tfhits[, bedCols], 
              paste0(filenameBase,'_TfNames.bed'), quote=F, row.names=F, col.names=F, sep='\t');
}



.checkGeneListEnrichment = function(dataset1, dataset2, refset, alt='two.sided') {
  data1name = as.character(match.call()[2]);
  data2name = as.character(match.call()[3]);
  incommon = sum(dataset1 %in% dataset2);
  outof = length(dataset2);
  innet = sum(dataset1 %in% refset);
  countdimnames = list(c("Y", "N", ""), c("Y", "N", ""));
  names(countdimnames)[[1]] = data2name;
  names(countdimnames)[[2]] = data1name;
  counts0 = matrix(ncol = 3, nrow = 3, dimnames = countdimnames);
  counts0[1,] = c(incommon, outof-incommon, outof);
  counts0[3,] = c(innet, length(refset)-innet, length(refset));
  counts0[2,] = counts0[3,] - counts0[1,];
  counts = counts0[-3,-3];
  test = fisher.test(counts, alternative = alt);
  out = list(counts0, test);
  return(out);
}

.checkGeneListEnrichmentList = function(dataset1, dataset2list, refset, thresh=.05, alt='two.sided', order=T) {			   	
  data1name = as.character(match.call()[2]);
  out=list(); 
  modpvals = as.data.frame(matrix(nrow=length(dataset2list), ncol=2, dimnames=list(names(dataset2list), c('pval', 'ratio'))));
  for (set in 1:length(dataset2list)) {
    if (suppressWarnings(all(sort(dataset1) == sort(dataset2list[[set]])))) { 
      next; 
    }
    temp = .checkGeneListEnrichment(dataset1, dataset2list[[set]], refset, alt=alt);
    names(dimnames(temp[[1]])) = c(names(dataset2list[set]), data1name);
    out[[set]] = temp;
    names(out)[set] = names(dataset2list[set]);
    modpvals[set, ] = c(signif(temp[[2]]$p.value, 3), signif(temp[[2]]$estimate, 3))
  }
  if(order){modpvals = modpvals[order(as.numeric(modpvals$pval)), ]};
  return(list(details = out, pvals = modpvals));
}

# assumes tfhitDf has names:
# seqnames, source, feature, start, end, absScore, relScore, strand, ID, TF, class, siteSeqs
# tfhitDf$seqnames is an interval with the format 'scaffold_0:100-200'
.computeTfHitAbsolutePositions = function (tfhitDf) {
  interval_split = strsplit(tfhitDf$seqnames[1], '[:-]')[[1]];
  interval_start = as.numeric(interval_split[2]);
  poscols = match(c('start','end'), names(tfhitDf));
  tfhitDf[, poscols] = tfhitDf[, poscols] + interval_start;
  tfhitDf$start = tfhitDf$start - 1;
  return(tfhitDf)
}

# assumes tfhitDf has names:
# seqnames, source, feature, start, end, absScore, relScore, strand, ID, TF, class, siteSeqs
.summarizeDmrTfHits = function (tfhitDf, stranded=F) {
  split_on = tfhitDf$TF;
  if (stranded) { split_on = paste0(split_on, tfhitDf$strand) }
  tfhitDf_summary = do.call('rbind', 
                            lapply(split(tfhitDf, split_on),
                                   function(x) data.frame(num=nrow(x),
                                                          mean.absScore=mean(x$absScore),
                                                          mean.relScore=mean(x$relScore))));
  tfhitDf_summary = as.data.frame(cbind(tfhitDf_summary, 
                                        num_x_absScore=tfhitDf_summary$num*tfhitDf_summary$mean.absScore));
  if (nrow(tfhitDf_summary)==1) {
    mean.rank = NA;
  } else {
    mean.rank = apply(apply(-tfhitDf_summary, 2, rank), 1, mean);
  }
  tfhitDf_summary = as.data.frame(cbind(tfhitDf_summary, mean.rank=mean.rank));
  return(tfhitDf_summary[order(tfhitDf_summary$mean.rank), ]);
}

.quantile.5 = function (vec) { return(quantile(vec, seq(.05,1,.05), na.rm=T)) }
.order2 = function(df, columns, decreasing=T) {
  if (is.character(columns)) {column=match(columns, colnames(df))}
  return(df[order(df[,columns], decreasing=decreasing), ])
}

.quantile.tails = function(vec) { return(quantile(vec, c(.01,.05,.1,.5,.9,.95,.99)))  }


# assumes both GRanges have mcols 
.findOverlapsNames = function (gr1, gr1nameCol, gr2=NULL, gr2nameCol=NULL, ...) {
  if (is.null(gr2)) {
    ovdf = as.data.frame(findOverlaps(gr1, ...));
  } else {
    ovdf = as.data.frame(suppressWarnings(findOverlaps(gr1, gr2, ...)));
  }
  ovdf$queryHits = mcols(gr1)[ovdf$queryHits, gr1nameCol];
  ovdf$subjectHits = mcols(gr2)[ovdf$subjectHits, gr2nameCol];
  return(ovdf);
}

# assumes gr1 has names
.findOverlapsWithSubjectMcols = function (gr1, gr2, ...) {
  ov = as.data.frame(suppressWarnings(findOverlaps(gr1, gr2, ...)));
  ov$queryHits = names(gr1)[ov$queryHits];
  ov = as.data.frame(cbind(queryHits=ov$queryHits,
                           as.data.frame(mcols(gr2[ov$subjectHits]))));
  return(ov)
}

.findOverlapsWithSelf = function (gr1, noRanges=F, ...) {
  ov0 = as.data.frame(suppressWarnings(findOverlaps(gr1, 
                                                    ignore.strand=T, drop.self=T, drop.redundant=T, 
                                                    ...)));
  ov = ov0;
  ov$queryHits = names(gr1)[ov$queryHits];
  ov$subjectHits = names(gr1)[ov$subjectHits];
  if (noRanges) {
    return(ov)
  } else {
    cat(paste0(nrow(ov), ' overlaps...\n'));
    ovdlist = list();
    mnum = ncol(mcols(gr1));
    for (i in 1:nrow(ov)) {
      if (i %% 100 == 0) { cat(i,'...') }
      ovdlist[[i]] = disjoin(gr1[unlist(ov0[i,])], ignore.strand=T)[2];
      inds = match(ov[i,], names(gr1));
      xr = as.data.frame(ranges(gr1[inds]));
      xm = as.data.frame(mcols(gr1[inds]));
      xs = as.vector(strand(gr1[inds]));
      tmc = cbind(xm[1,], xr[1,1:3], xs[1], xm[2,], xr[2,1:3], xs[2]);
      names(tmc)[c((mnum+4),ncol(tmc))] = rep('strand',2);
      names(tmc) = c(paste0(names(tmc)[1:(mnum+4)],'.1'),
                     paste0(names(tmc)[(mnum+5):(ncol(tmc))],'.2'));
      mcols(ovdlist[[i]]) <- tmc;
    }
    return(list(ov_raw=ov0,ov=ov,ov_ranges=ovdlist));
  }
}

# assumes gene name is in first column of mcols for gr1
.getOVgenes = function(gr1,gr2) {
  return(suppressWarnings(mcols(gr1)[unique(as.data.frame(findOverlaps(gr1,gr2))$queryHits), 1]));
}


.verboseBoxplotColumns = function (mat, factor, mfrow=NULL, hline=NULL, ...) {
  if (is.null(mfrow)) {
    mfrow = c(2, ceiling(ncol(mat)/2));
  }
  par(mfrow=mfrow);
  for (i in 1:ncol(mat)) {
    WGCNA::verboseBoxplot(mat[,i], as.factor(factor), ylab=colnames(mat)[i], xlab='', ...);
    if (is.numeric(hline)) {
      abline(h=hline);
    }
  }
}

# assumes gr has names
.distance.named = function (gr, gr_interval) {
  d = distance(gr, gr_interval);
  names(d) = names(gr);
  return(d[!is.na(d)]);
}

.getCGforDNAString = function (DNAString) {
  freqs = alphabetFrequency(DNAString);
  GC = sum(freqs[match(c('C','G'), names(freqs))]) / sum(freqs);
  return(GC);
}

.generateWindows = function (start, end, size, step) {
  #startvec = seq(start, end - size, step);
  startvec = seq(start, end, step);
  endvec = startvec + size - 1;
  endvec[endvec > end] = end;
  return(cbind(start=startvec, end=endvec))
}

.generateDNAStringSubstrings = function (DNAString, window_size, step_size) {
  seqlen = length(DNAString)
  starts = seq(1, seqlen, step_size); 
  intervals = cbind(start=starts, end=(starts+window_size-1)); 
  intervals[intervals[,2] > seqlen, 2] = seqlen; 
  subs = apply(intervals, 1, 
               function(f) subseq(DNAString, start=f[1], end=f[2]));
  return(list(windows=intervals, subseqs=subs));
}

.getGCcontentForDNAStringWindows = function (DNAString, window_size, step_size) {
  x = .generateDNAStringSubstrings(DNAString, window_size, step_size);
  s = sapply(x$subseqs, .getCGforDNAString);
  return(list(GC=s, subseqs=x));
}

.h = function (df, n=10) { return(df[1:n,1:n]) }

.screenGOFunctionAllOntologies = function (idlist, bgvec, fdrmethvec=c('BY','BH'), fdrthvec=c(.01, .05, .1)) {
  golist = list();
  for (fdrmeth in fdrmethvec) {
    cat(paste0('------------------------- ', fdrmeth,  ' -------------------------'));
    for (fdrth in fdrthvec) {
      cat(paste0('------------------------- ', fdrth,  ' -------------------------'));
      for (i in 1:length(idlist)) {
        cat(paste0('------------------------- ', names(idlist)[i],  ' -------------------------'));
        thisres = .GOFunctionAllOntologies(idlist[[i]], bgvec, fdrmethod=fdrmeth, fdrth=fdrth);
        thisname = paste0(fdrmeth,'.', gsub('0.','',fdrth,fixed=T));
        if (nrow(thisres) > 0) {
          golist[[names(idlist)[i]]][[thisname]] = .addGenesToGOFunctionResults(thisres, modulegenes=idlist[[i]]);
        } else {
          golist[[names(idlist)[i]]][[thisname]] = thisres;
        }
      }
    }
  }
  return(golist);
}

.buildEntrezListsFromAbGeneLists = function (abgenelist, abHsMap, limitnum, ...) {
  entrezlists = list();
  for (f in 1:length(abgenelist)) {
    tmap = .mapAbGenesHsEntrez(unique(unlist(abgenelist[[f]])), abHsMap, 'ab');
    entrezlists[[names(abgenelist)[f]]] = .limitVectorsInList(tmap, num=limitnum, smartUnlist=T, ...)$newvec;
  }
  return(entrezlists);
}

.getGenesByTerm = function (GOFunctionOutput, geneCol, geneDelim, termCol, originalAbGenes=NULL) {
  x = strsplit(GOFunctionOutput[,geneCol], geneDelim, fixed=T);
  names(x) = GOFunctionOutput[,termCol];
  if (is.vector(originalAbGenes)) {
    x = lapply(x, function(f) originalAbGenes[originalAbGenes %in% f]);
  }
  return(x)
}

.getTermsByGene = function (getGenesByTermOutput) {
  x = names(unlist(getGenesByTermOutput));
  x = strsplit(x, '.', fixed=T);
  out = list();
  for (i in x) {
    if (i[2] %in% names(out)) {
      out[[i[2]]] = c(out[[i[2]]], i[1]);
    } else {
      out[[i[2]]] = i[1];
    }
  }
  return(out)
}

.parseGOhitGenes = function (GOFunctionOutput, geneCol, geneDelim, termCol, originalAbGenes=NULL) {
  genesByTerm = .getGenesByTerm(GOFunctionOutput, geneCol, geneDelim, termCol, originalAbGenes=originalAbGenes);
  termsByGene = .getTermsByGene(genesByTerm);
  return(list(byTerm=genesByTerm, byGene=termsByGene));
}

.buildGeneTermTable = function (parseGOhitGenesOutput) {
  bygene = parseGOhitGenesOutput$byGene;
  byterm = parseGOhitGenesOutput$byTerm;
  df = as.data.frame(matrix(nrow=length(bygene), ncol=length(byterm), 
                            dimnames=list(names(bygene),names(byterm))));
  for (g in 1:length(bygene)) {
    df[g, ] = as.numeric(names(df) %in% bygene[[g]]);
  }
  return(df)
}

.getDmrsAndAbGenesForTermEntrezGenes = function (GOFunctionAllOntologiesOutput, inputGenesToGOFunction, dmrsOVgene,
                                                 geneCol=8, geneSep=',', termCol=3) {
  go = GOFunctionAllOntologiesOutput;
  inputgenes = inputGenesToGOFunction;
  go.parsed = .parseGOhitGenes(thisGO, geneCol, geneSep, termCol, inputgenes);  
  abgenesByTerm = lapply(go.parsed$byTerm, 
                         function(f) sapply(strsplit(names(f),'__'), 
                                            function(ff) ff[1]));
  dmrsOVgeneByTerm = lapply(abgenesByTerm, 
                            function(f) dmrsOVgene[match(f, dmrsOVgene$mcols.geneSym), ]);
  for (tt in 1:length(dmrsOVgeneByTerm)) {
    dmrsOVgeneByTerm[[tt]] = as.data.frame(cbind(dmrsOVgeneByTerm[[tt]],
                                                 entrez=go.parsed$byTerm[[tt]]));
  }
  return(dmrsOVgeneByTerm);
}




.screen_enrichGOAllOntologies = function (idlist, bgvec, OrgDb=org.Hs.eg.db, fdrmethvec=c('BY','BH'), fdrthvec=c(.01, .05)) {
  golist = list();
  for (fdrmeth in fdrmethvec) {
    cat(paste0('------------------------- ', fdrmeth,  ' -------------------------'));
    for (fdrth in fdrthvec) {
      cat(paste0('------------------------- ', fdrth,  ' -------------------------'));
      for (i in 1:length(idlist)) {
        cat(paste0('------------------------- ', names(idlist)[i],  ' -------------------------'));
        thisres = .enrichGOAllOntologies(idlist[[i]], OrgDb=OrgDb, refGenes=bgvec, pAdjustMethod=fdrmeth, pvalueCutoff=fdrth);
        #thisres = .GOFunctionAllOntologies(idlist[[i]], bgvec, fdrmethod=fdrmeth, fdrth=fdrth);
        thisname = paste0(fdrmeth,'.', gsub('0.','',fdrth,fixed=T));
        # if (nrow(thisres) > 0) {
        #   golist[[names(idlist)[i]]][[thisname]] = .addGenesToGOFunctionResults(thisres, modulegenes=idlist[[i]]);
        # } else {
        #   golist[[names(idlist)[i]]][[thisname]] = thisres;
        # }
        golist[[names(idlist)[i]]][[thisname]] = thisres;
      }
    }
  }
  return(golist);
}