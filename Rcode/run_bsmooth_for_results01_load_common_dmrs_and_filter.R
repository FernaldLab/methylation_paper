# clear workspace, set options and working directory, load functions
rm(list=ls()); options(stringsAsFactors=F); setwd('~/Documents/_BS-seq_analysis/');
library('GenomicRanges');
source('/Volumes/fishstudies-1/_code/run_bsmooth_functions.R');

# # load common dmrs
newprefix = 'CpGcombined.CG_4xCov_ns70_50_25hAll';
pattern = 'mg1e+08_fits_4xCov_Ttests_DMRs_cut0.95_mg300_tstat.corrected_n3_meanDiff0.1.RData';
# filebase = gsub('.RData', '', paste0(newprefix, pattern), fixed=T);
# load(paste0(filebase, 'exactAnd100bpOverlaps.RData'));
# 
# # combine overlapping dmrs and exact matches
# dmrs0 = as.data.frame(rbind(ovdmrs, mdmrs));
load('WORKSPACE_CpGcombined.CG_4xCov_ns70_50_25hAllmg1e+08_fits_4xCov_Ttests_DMRs_cut0.95_mg300_tstat.corrected_n3_meanDiff0.1.RData');
dmrs0 = DMRs25_50_70; rm(list=grep('70$',ls(),value=TRUE));
dmrs0ints = paste0(dmrs0$chr, ':', dmrs0$start, '-', dmrs0$end);
names(dmrs0ints) = dmrs0$direction;
.write_bedfile_for_intervals(dmrs0ints, 'dmrs0_1872_DMRs25_50_70.bed', returnDataFrame=F);

# load dataframe containing hand curation of dmrs
hc2 = read.table(gsub('RData','txt',paste0('curatedALL_',newprefix,pattern)), header=T, sep='\t');
hc2GR = GRanges(seqnames=hc2$chr, ranges=IRanges(start=hc2$start, end=hc2$end), mcols=hc2$status);
hc2ints = paste0(hc2$chr, ':', hc2$start, '-', hc2$end);
if (any(hc2ints[match(dmrs0ints,hc2ints)] != dmrs0ints)) { stop('hc2 and dmrs0 don\'t match') }

# get ids of dmrs that passed curation and filter
goodstatus = c('ok','ok_s','ok_d');
goodhc2ints = hc2ints[hc2$status %in% goodstatus];

dmrs = dmrs0[dmrs0ints %in% goodhc2ints, ];
dmrsints = paste0(dmrs$chr, ':', dmrs$start, '-', dmrs$end);
rownames(dmrs) = dmrsints;
names(dmrsints) = dmrs$direction;
.write_bedfile_for_intervals(dmrsints, 'dmrs_1034_hc.ok_ok_s_ok_d.bed', returnDataFrame=F);

mcolCols = which(names(dmrs)=='n'):ncol(dmrs);
dmrsGR = GRanges(seqnames=dmrs$chr, ranges=IRanges(start=dmrs$start, end=dmrs$end), mcols=dmrs[,mcolCols]);
names(dmrsGR) = dmrsints;

rm(list=grep('dmrs',ls(),invert=T,value=T));
save(list=ls(), file='run_bsmooth_for_results01_load_common_dmrs_hc_filterWORKSPACE.RData');

# load null dmrs
load('WORKSPACE_run_bsmooth_generate_nulls.RData');

# find and filter out dmrs that are overlapped more than 50% by a null dmr
allnullints = unique(unlist(nullints));
.write_bedfile_for_intervals(allnullints, 'allnullints.bed', returnDataFrame=F);

# allnullintsGR = .buildGR_for_intervals(allnullints);
# nov = list()
# for (d in 1:nrow(dmrs)) {
#   thisw = dmrs$width[d];
#   nov[[dmrsints[d]]] = as.data.frame(findOverlaps(dmrsGR[d], allnullintsGR, minoverlap=(thisw/2)))
# }; rm(d,thisw);
# 
# dmrs = dmrs[!(dmrsints %in% names(nov)[sapply(nov, nrow)!=0]), ];

dmrs_filtered_null.9 = .filter_nullints_fromDMRs(dmrs, allnullints, .9);
dmrs_filtered_null.5 = .filter_nullints_fromDMRs(dmrs, allnullints, .5);
dmrs_filtered_null.33 = .filter_nullints_fromDMRs(dmrs, allnullints, .33);
dmrs_filtered_null.1 = .filter_nullints_fromDMRs(dmrs, allnullints, .1);

towrite = rownames(dmrs_filtered_null.9);
names(towrite) = dmrs_filtered_null.9$direction;
.write_bedfile_for_intervals(towrite, 'dmrs_769_hc.ok_ok_s_ok_d_null.9.bed', returnDataFrame=F);

towrite = rownames(dmrs_filtered_null.5);
names(towrite) = dmrs_filtered_null.5$direction;
.write_bedfile_for_intervals(towrite, 'dmrs_709_hc.ok_ok_s_ok_d_null.5.bed', returnDataFrame=F);

towrite = rownames(dmrs_filtered_null.33);
names(towrite) = dmrs_filtered_null.33$direction;
.write_bedfile_for_intervals(towrite, 'dmrs_673_hc.ok_ok_s_ok_d_null.33.bed', returnDataFrame=F);

towrite = rownames(dmrs_filtered_null.1);
names(towrite) = dmrs_filtered_null.1$direction;
.write_bedfile_for_intervals(towrite, 'dmrs_604_hc.ok_ok_s_ok_d_null.1.bed', returnDataFrame=F);

rm(list=grep('dmrs_filtered',ls(),value=T,invert=T));
save(list=ls(), file='run_bsmooth_for_results01_load_common_dmrs_hc_and_null_filterWORKSPACE.RData');
