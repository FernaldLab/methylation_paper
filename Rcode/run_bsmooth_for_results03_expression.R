rm(list=ls()); options(stringsAsFactors=F); 
setwd('~/Documents/_BS-seq_analysis/');
library(DESeq2);
load('annotations_and_scaffold_stats_after03_get_overlaps.RData');
load('GRanges_objects_after03_get_overlaps.RData');
load('dmrs_and_overlaps_after03_get_overlaps.RData');

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# load gene expression data
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# load estimated count data from kallisto
exprdir = '/Volumes/fishstudies-1/_LYNLEY_RNAseq/data/_kallisto_092915/';
nd1exp = read.table(paste0(exprdir, 'ATCACG/abundance.tsv'), header=T, sep='\t');
rownames(nd1exp) = sapply(strsplit(nd1exp$target_id, '|', fixed=T), function(f) f[length(f)]);
nd2exp = read.table(paste0(exprdir, 'TGACCA/abundance.tsv'), header=T, sep='\t');
rownames(nd2exp) = sapply(strsplit(nd2exp$target_id, '|', fixed=T), function(f) f[length(f)]);
d1exp = read.table(paste0(exprdir, 'CGATGT/abundance.tsv'), header=T, sep='\t');
rownames(d1exp) = sapply(strsplit(d1exp$target_id, '|', fixed=T), function(f) f[length(f)]);
d2exp = read.table(paste0(exprdir, 'TTAGGC/abundance.tsv'), header=T, sep='\t');
rownames(d2exp) = sapply(strsplit(d2exp$target_id, '|', fixed=T), function(f) f[length(f)]);

# make count matrix, columns are subjects, rows are transcripts
counts0 = as.data.frame(cbind('3157_TENNISON'=nd1exp$est_counts, 
                              '3165_BRISCOE'=d1exp$est_counts, 
                              '3581_LYNLEY'=d2exp$est_counts, 
                              '3677_MONK'=nd2exp$est_counts));

an = annotations_and_scaffold_stats$an;

# sum across transcripts to get gene-level expression values
rownames(counts0) = rownames(nd1exp);
counts0 = counts0[match(an$lookup$transcript_id,rownames(counts0)), ];
counts0l = split(counts0,an$lookup$gene);
counts = as.data.frame(t(sapply(counts0l, function(f) apply(f, 2, sum))));
rm(counts0, counts0l); gc();


# ----------------------------------------------------------------------------------------------------------
# compute differential expression and save results
# ----------------------------------------------------------------------------------------------------------
# make design matrix
coldata = as.data.frame(matrix(c('ND','D','D','ND'), ncol=1));
rownames(coldata) = names(counts);
names(coldata) = 'group';

# run DESeq2
dds = DESeqDataSetFromMatrix(countData=round(counts), colData=coldata, design = ~ group); 
dds = DESeq(dds, parallel=TRUE); 
res_raw = results(dds, contrast=c('group','D','ND'));
res = as.data.frame(res_raw);
res = as.data.frame(cbind(res[,c(1,2,6)], 
                          an$gffGenesDF[match(rownames(res), an$gffGenesDF$geneSym), ],
                          brGC=mcols(grs$brGR)$CG[match(rownames(res), names(grs$brGR))]));
rm(first_time); # don't know what this is but gets created in workspace by DESeq functions
####################

dmrs = dmrs_and_overlaps$dmrs;

# subset deseq results to expressed genes, DE genes, and genes on dmr-containing scaffolds
resexpr = subset(res, baseMean>0);
res.1 = subset(res, padj<.1);
resexprDmrScaffolds = subset(resexpr, seqnames %in% unique(dmrs$chr));

# make res version with genes from all scaffolds that were smoothed
# save
expr_raw = list(nd1=nd1exp, nd2=nd2exp, d1=d1exp, d2=d2exp);
deseq_input = list(counts=counts, data=dds, coldata=coldata);
deseq_results = list(raw=res_raw, raw_df=res, DE=res.1, expressed=resexpr, expressed_dmrScaffolds=resexprDmrScaffolds);

save(expr_raw, deseq_input, deseq_results, file='deseq_methylationDvsND_res.RData');
####################

# add dmr information
dmrsOV = dmrs_and_overlaps$dmrsOV;
dmrsOVsimple = dmrs_and_overlaps$dmrsOVsimple;
dmrsOVbyGene = dmrs_and_overlaps$dmrsOVbyGene;
alldmrgenes = unique(unlist(dmrsOVsimple));

res = deseq_results$raw_df;

res = .addDmrInfoToExpression(res, dmrsOVbyGene);

with_te = annotations_and_scaffold_stats$teOVwith_features;
for (i in 1:length(with_te)) {
  res = as.data.frame(cbind(res, 
                            rownames(res) %in% with_te[[i]]));
  names(res)[ncol(res)] = paste0(names(with_te)[i],'.teOV');
}; rm(i);

res.1 = subset(res, padj<.1);
resexpr = subset(res, baseMean>0);
resexprSmoothed = subset(resexpr, 
                         seqnames %in% rownames(annotations_and_scaffold_stats$slGeneStats_in70))
resexprSmoothedDmrs = subset(resexprSmoothed, seqnames %in% unique(dmrs$chr));

save(res, res.1, resexpr, resexprSmoothed, resexprSmoothedDmrs, 
     file='deseq_methylationDvsND_res_withDmrInfo_noDESeq_input_or_raw_expression.RData');

#####################






















# get expression for genes with dmrs in introns
#thesegenes = names(dmrsOVbyGene$intron)

#thisOVbyfeat = thisOVbyfeat[sapply(thisOVbyfeat, function(f) all(f$direction=='hypo'))]
thismap = genebrmaps;
thisresexpr = resexpr_in70dmrsOV;

# test intron hits
thismap = subset(thismap, mcols.geneSym %in% names(dmrsOVbyGene$intron));
thisfg = fg[match(thismap$queryHits, names(fg))];
thisres = thisresexpr[match(thismap$mcols.geneSym, rownames(thisresexpr)), ]
rownames(thisres) = thismap$mcols.geneSym

verboseBoxplot(thisres$baseMean, 
               sapply(strsplit(rownames(genebrmaps1),'.',fixed=T), 
                      function(f) f[1]), 
               ylab='', ylim=c(0,5000))



x = thisfg %in% names(table(thisfg))[table(thisfg) > 9];
thisfg = thisfg[x]
thisres = thisres[x, ]
thismap = thismap[x,]

##
xmaps = subset(genebrmaps, !grepl(':',fg))
thesegenes = unique(genebrmaps$mcols.geneSym)
thisres = thisresexpr[rownames(thisresexpr) %in% thesegenes,];
thesedmrs = thisOVbyfeat[match(rownames(thisres), thesegenes)];

thesefg = lapply(lapply(thesedmrs, rownames), 
                 function(f) unique(fg[match(f, names(fg))]));

singlefg = sapply(thesefg,length)==1;
thesefg = unlist(thesefg[singlefg]);
thisres = thisres[singlefg, ];

x = thesefg %in% c('intron', 'cds:exon:intron', 'intron:te', 'intron:z', 'intron:br.diff', 'intron:br3.diff')
thesefg = thesefg[x];
thisres = thisres[x, ];

#verboseBoxplot(thisres$baseMean, sapply(thesefg, function(f) all(f=='intron')))
verboseBoxplot(thisres$baseMean, thesefg=='intron')
