# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# load burtoni genome information and features
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# load a list with the following elements:
#  $gffGenes - gene gff dataframe
#  $gffGenesGR - gene gff GRanges
#  $gff - full gff in dataframe
#  $recipBL - dataframe mapping burtoni genes to other fish and human ids
#  $gois - vector of ids for typical Fernald lab genes of interest
#  $lookup - dataframe mapping burtoni transcript and protein ids to genes
an = .loadBurtoniGenesApr2016();

# load a list containing other types of genomic features
an2 = .loadBurtoniOtherFeatures('/Volumes/fishstudies-1/_Burtoni_annotations/');

# load a dataframe with mappings of burtoni genes to human ids that aren't in an$recipBL
handannos = read.table('~/Documents/_BS-seq_analysis/dmr_handAnnotations_hsaEntrez.txt', sep='\t', header=F, quote='');
handannos = handannos[!duplicated(handannos[,2]), ];
handannos = handannos[grep('^[0-9]', handannos[,1]), ];
names(handannos) = c('hsaHomologEntrez','gene','description');

# combine reciprocal BLAST results and hand annotations for mapping burtoni genes to human entrez ids
annoCombo = .combineHandAndMachineAnnos(an$recipBL, handannos);
abHsMap = strsplit(annoCombo$new$hsaHomologEntrez, ',');
names(abHsMap) = annoCombo$new$gene;
abHsMap = abHsMap[!sapply(abHsMap, function(f) all(is.na(f)))];

# load scaffold lengths into a named vector
sl = read.table('/Volumes/fishstudies-1/_Burtoni_genome_files/scaffold_lengths'); 
sl0 = sl; sl = sl0[, 2]; names(sl) = sl0[, 1];

# build table of statistics about genes on each scaffold
slGeneStats = as.data.frame(matrix(nrow=length(sl), ncol=7));
statnames = c('length','numGenes','avgGeneWidth','avgGeneGC','avgNumIsoforms','pctBpGene','genesPerMb');
dimnames(slGeneStats) = list(names(sl), statnames);
slGeneStats$length = sl;

# count genes on each scaffold
tmpgenenum = table(an$gffGenesDF$seqnames);
tmprows = match(rownames(slGeneStats), names(tmpgenenum));
slGeneStats$numGenes[which(!is.na(tmprows))] = tmpgenenum[na.omit(tmprows)];
rm(tmpgenenum, tmprows);
slGeneStats$numGenes[which(is.na(slGeneStats$numGenes))] = 0;

# get other statistics about genes by scaffold
for (i in which(slGeneStats$numGenes != 0)) {
  sc = rownames(slGeneStats)[i];
  if (i %% 1000 == 0) {cat(sc,'...')}
  scrows = an$gffGenesDF$seqnames == sc;
  slGeneStats$avgGeneWidth[i] = mean(an$gffGenesDF$width[scrows]);
  slGeneStats$avgGeneGC[i] = signif(mean(an$gffGenesDF$GC[scrows]), 4);
  slGeneStats$avgNumIsoforms[i] = signif(mean(an$gffGenesDF$isoforms[scrows]), 4);
  slGeneStats$pctBpGene[i] = signif(slGeneStats$avgGeneWidth[i] / slGeneStats$length[i], 4);
  slGeneStats$genesPerMb[i] = slGeneStats$numGenes[i] / slGeneStats$length[i] * 1000000
}; rm(i,sc,scrows);

#
.verboseScatterplotAllColumnPairs(slGeneStats, mfrow=c(3,7), abline.col='red',
                                  bg='grey', col='black', pch=21, cex=1.8);

#
.verboseScatterplotAllColumnPairs(an$gffGenesDF[,c(4,9,10)], mfrow=c(1,3), abline.col='red',
                                  bg=WGCNA::labels2colors(an$gffGenesDF$strand), col='black', 
                                  pch=21, cex=1.8);
.verboseScatterplotAllColumnPairs(an$gffGenesDF[,c(4,9,10)], mfrow=c(1,3), abline.col='red',
                                  bg=WGCNA::labels2colors(an$gffGenesDF$strand), col='black', 
                                  pch=21, cex=1.8,
                                  corOptions='method="s"');


# load dataframe of non-coding regions conserved from drerio
zcne = read.table('~/Documents/_annotationsDec2015/hBurtoni_zCNE_alignments.out.bed', sep='\t');




# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# make GRanges objects
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------


# gene bodies, transcription start sites, basal regulatory regions
geneGR = an$gffGenesGR;
# tssGR = .getFeatureStartsFromGR(geneGR, stranded=T)$gr;
basalReg = .getWindowsAroundGeneStarts(gr=geneGR, up=5000, down=1000, sortOutput=F);
basalReg = .fixInvalidPositions(basalReg$gr, chrLengthVec=sl);
brGR = basalReg$gr; rm(basalReg);

# 3' basal regulatory region
br3primeGR = .getWindowsAroundGeneEnds(gr=geneGR, down=5000, up=1000, sortOutput=F);
br3primeGR = .fixInvalidPositions(br3primeGR$gr, chrLengthVec=sl)$gr;

# make exon GR
gffExon0 = subset(an$gff, V3=='exon');
#gffExon = gffExon0[!duplicated(paste0(gffExon0$V1, gffExon0$V4, gffExon0$V5)), ];
gffExon = gffExon0
exonGR = GRanges(seqnames=gffExon$V1, 
                 ranges=IRanges(start=gffExon$V4, end=gffExon$V5), 
                 strand=gffExon$V7, 
                 mcols=as.data.frame(cbind(.parseGffMetaCol(gffExon, pattern='gene='),
                                           gffExon$V11,
                                           gffExon$V18)));
rm(gffExon, gffExon0); gc();

# make cds GR
gffCDS0 = subset(an$gff, V3=='CDS');
#gffCDS = gffCDS0[!duplicated(paste0(gffCDS0$V1, gffCDS0$V4, gffCDS0$V5)), ];
gffCDS = gffCDS0;
cdsGR = GRanges(seqnames=gffCDS$V1, 
                ranges=IRanges(start=gffCDS$V4, end=gffCDS$V5), 
                strand=gffCDS$V7, 
                mcols=as.data.frame(cbind(.parseGffMetaCol(gffCDS, pattern='gene='),
                                          gffCDS$V11,
                                          gffCDS$V18)));
rm(gffCDS, gffCDS0); gc();

# make mRNA GR 
gffmrna0 = subset(an$gff, V3=='mRNA');
#gffmrna = gffmrna0[!duplicated(paste0(gffmrna0$V1, gffmrna0$V4, gffmrna0$V5)), ];
gffmrna = gffmrna0;
mrnaGR = GRanges(seqnames=gffmrna$V1, 
                 ranges=IRanges(start=gffmrna$V4, end=gffmrna$V5), 
                 strand=gffmrna$V7, 
                 mcols=as.data.frame(cbind(.parseGffMetaCol(gffmrna, pattern='gene='),
                                           gffmrna$V11,
                                           gffmrna$V18)));
rm(gffmrna0, gffmrna); gc();

# make intron GR
load('~/Documents/_annotationsDec2015/intron.gff.RData');
x.int_rnanum = .parseGffMetaCol(x.int, pattern='Parent=');
gffIntron = as.data.frame(cbind(x.int, gene=an$lookup$gene[match(x.int_rnanum, an$lookup$rna_num)]));
intronGR = GRanges(seqnames=gffIntron$V1, 
                   ranges=IRanges(start=gffIntron$V4, end=gffIntron$V5), 
                   strand=gffIntron$V7, 
                   mcols=as.data.frame(cbind(geneSym=gffIntron$gene,
                                             gffIntron$V9)));
rm(x.int, x.int_rnanum, gffIntron); gc();

# transposons, microRNAs, conserved non-coding regions
teGR = an2$gr$te;
miGR = an2$gr$mi;
zranges = IRanges(start=zcne[,2], end=zcne[,3]);
zGR = GRanges(seqnames=zcne[,1], ranges=zranges, strand=zcne[,6], mcols=zcne[,4:5]);
rm(zranges)

save(handannos, slGeneStats, sl, sl0, abHsMap, an, an2, annoCombo, zcne,
     br3primeGR, brGR, cdsGR, exonGR, geneGR, intronGR, miGR, mrnaGR, teGR, zGR,
     file='~/Documents/_annotationsDec2015/WORKSPACE_annotations_prep_for_dmrs_wgcna.RData');








######################
#### functions

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

#################

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






.parseGffMetaCol = function(gff, metaCol=9, delim=';', pattern='Dbxref=GeneID:') {
  metaSplit = strsplit(gff[, metaCol], delim);
  check = sapply(metaSplit, function(f) any(grepl(pattern, f)))
  if (!all(check) | sum(check)!=nrow(gff)) {
    stop('all lines in gff do not contain arg pattern') 
  }
  res = grep(pattern, unlist(metaSplit), value=TRUE);
  return(gsub(pattern, '', res));
}
