# -----------------------------------------------------------------------------------------------
# setup
# -----------------------------------------------------------------------------------------------
rm(list=ls()); options(stringsAsFactors=F); 
library(org.Hs.eg.db);
library(GOFunction);
library(clusterProfiler);
library(VennDiagram);

# get list with burtoni-human mappings (abHsMap)
load('abHsMap.RData');

# get list with vectors of burtoni genes (abgenelist)
load('example_gene_lists.RData');

# build background for enrichment tests
bgGO = unique(unlist(abHsMap));

# make list for mapping GOFunction terms to genes, need for .getGenesForGOTerms()
org.Hs.egGO2ALLEGSlist = as.list(org.Hs.egGO2ALLEGS);

# convert burtoni genes to human entrez ids
entrez_lists = .buildEntrezListsFromAbGeneLists(abgenelist, abHsMap, limitnum=1, method='first');

# -----------------------------------------------------------------------------------------------
# run GO functions
# -----------------------------------------------------------------------------------------------
fdrmethvec = c('BH','BY');
fdrthvec = c(.01,.05,.1);

# -----------------------------------------------------------------------------------------------
# use GOFunction for more conservative testing
go1.0 = .screenGOFunctionAllOntologies(entrez_lists, bgGO, fdrmethvec=fdrmethvec, fdrthvec=fdrthvec);

# parse some results from GOFunction and cluster GO terms with common genes
go1 = go1.0$d$BY.1;
go1parsed  = .parseGOhitGenes(go1, 8, ',', 3, entrez_lists$d);
go1parsedDf = .buildGeneTermTable(go1parsed);
names(go1parsedDf) = paste0(1:ncol(go1parsedDf), '-', names(go1parsedDf));
heatmap(1-cor(go1parsedDf), symm=T);

# -----------------------------------------------------------------------------------------------
# use enrichGO for less conservative testing
go2.0 = .screen_enrichGOAllOntologies(entrez_lists, bgGO, fdrmethvec=fdrmethvec, fdrthvec=fdrthvec);

# parse some results from enrichGO and cluster terms
go2 = go2.0$d$BY.1;
go2parsed  = .parseGOhitGenes(go2, 9, '/', 3, entrez_lists$d);
go2parsedDf = .buildGeneTermTable(go2parsed);
names(go2parsedDf) = paste0(1:ncol(go2parsedDf), '-', names(go2parsedDf));
heatmap(1-cor(go2parsedDf), symm=T);

# -----------------------------------------------------------------------------------------------
# make venn diagrams to compare results, also a function defined for triple venns
ov_2 = .draw_pairwise_venn_from_groups(go2.0$d$BY.1$Description, go2.0$nd$BY.1$Description, 
									   groupNames=c('d','nd'), 
									   fill=c('lightblue','greenyellow'), col='grey');
									   
ov_4 = .draw_quad_venn_from_groups(go2.0$d$BY.1$ID, go2.0$nd$BY.1$ID, go1.0$d$BY.1$goid, go1.0$nd$BY.1$goid,
								   fill=c('lightblue','blue','greenyellow','green'), col='grey');

# -----------------------------------------------------------------------------------------------
# kegg pathways
kegg1 = .screen_enrichKEGG(entrez_lists, bgGO, fdrmethvec=c('BH','BY'), fdrthvec=c(.01,.05,.1));

# ===============================================================================================
# functions
# ===============================================================================================

# -----------------------------------------------------------------------------------------------
# converting burtoni ids to human and gene ontology
# -----------------------------------------------------------------------------------------------
.buildEntrezListsFromAbGeneLists = function (abgenelist, abHsMap, limitnum, ...) {
  entrezlists = list();
  for (f in 1:length(abgenelist)) {
    tmap = .mapAbGenesHsEntrez(unique(unlist(abgenelist[[f]])), abHsMap, 'ab');
    entrezlists[[names(abgenelist)[f]]] = .limitVectorsInList(tmap, num=limitnum, smartUnlist=T, ...)$newvec;
  }
  return(entrezlists);
}

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

.screenGOFunctionAllOntologies = function (idlist, bgvec, fdrmethvec=c('BY','BH'), fdrthvec=c(.01, .05, .1)) {
  golist = list();
  for (fdrmeth in fdrmethvec) {
    cat(paste0('------------------------- ', fdrmeth,  ' -------------------------\n'));
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

.addGenesToGOFunctionResults = function (GOFunctionAllOntologiesOutput, maps=org.Hs.egGO2ALLEGSlist, modulegenes) {
  genelist = .getGenesForGOTerms(GOFunctionAllOntologiesOutput$goid, genes=modulegenes);
  new = as.data.frame(cbind(GOFunctionAllOntologiesOutput, genes=rep(NA, nrow(GOFunctionAllOntologiesOutput))));
  for (row in 1:nrow(new)) {
    tgenes = genelist[[match(new$goid[row], names(genelist))]];
    new$genes[row] = paste(tgenes, collapse=',');
  }
  return(new);
}

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

.screen_enrichGOAllOntologies = function (idlist, bgvec, OrgDb=org.Hs.eg.db, fdrmethvec=c('BY','BH'), fdrthvec=c(.01, .05), ...) {
  golist = list();
  for (fdrmeth in fdrmethvec) {
    cat(paste0('------------------------- ', fdrmeth,  ' -------------------------\n'));
    for (fdrth in fdrthvec) {
      cat(paste0('------------------------- ', fdrth,  ' -------------------------'));
      for (i in 1:length(idlist)) {
        cat(paste0('------------------------- ', names(idlist)[i],  ' -------------------------'));
        thisres = .enrichGOAllOntologies(idlist[[i]], OrgDb=OrgDb, refGenes=bgvec, pAdjustMethod=fdrmeth, pvalueCutoff=fdrth, ...);
        thisname = paste0(fdrmeth,'.', gsub('0.','',fdrth,fixed=T));
        golist[[names(idlist)[i]]][[thisname]] = thisres;
      }
    }
  }
  return(golist);
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

.screen_enrichKEGG = function (idlist, bgvec, fdrmethvec=c('BY','BH','fdr'), fdrthvec=c(.01, .05, .1), ...) {
	kegglist = list();
	for (fdrmeth in fdrmethvec) {
    cat(paste0('\n------------------------- ', fdrmeth,  ' -------------------------'));
    for (fdrth in fdrthvec) {
      cat(paste0('\n  ', fdrth,  '...\n'));
      for (i in 1:length(idlist)) {
        cat(paste0('     ', names(idlist)[i], '...'));
        thisres = as.data.frame(enrichKEGG(gene=idlist[[i]], universe=bgvec, pAdjustMethod=fdrmeth, pvalueCutoff=fdrth, ...));
        thisname = paste0(fdrmeth,'.', gsub('0.','',fdrth,fixed=T));
        kegglist[[names(idlist)[i]]][[thisname]] = thisres;
      }
    }
  }
  return(kegglist);
}

.parseGOhitGenes = function (GOFunctionOutput, geneCol, geneDelim, termCol, originalAbGenes=NULL) {
  genesByTerm = .getGenesByTerm(GOFunctionOutput, geneCol, geneDelim, termCol, originalAbGenes=originalAbGenes);
  termsByGene = .getTermsByGene(genesByTerm);
  return(list(byTerm=genesByTerm, byGene=termsByGene));
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

# -----------------------------------------------------------------------------------------------
# compare lists of genes, terms, etc.
# -----------------------------------------------------------------------------------------------
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
                 ...);
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
                      paste0(groupNames[1:3],collapse=':'));
  return(list(overlaps=overlaps, uniques=uniques, groupNames=groupNames));
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
                 ...);
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
                      paste0(groupNames[1:4],collapse=':'));
  return(list(overlaps=overlaps, uniques=uniques, groupNames=groupNames));
}

