#' @title DiffMetaScreener: a multi-level meta-analysis framework robustly screens directional signatures based on difference analyses.
#' @param expression_accession_vector A vector of all the sub-dataset names.
#' @param comparisons A file name is needed for each sub-dataset that distinguishes between high(activation) and low(inhibition) biological activity. File column 1: Accession, column 2: Group1, column 3: Group2, column 4: Group1_Status, column 5: Group2_Status.
#' @param geneSets_gmt Genesets to be screened.(Please make sure that every gene set in there is unique, otherwise there will be severe covariance.)
#' @param enrichment_method Selection of enrichment methods, multiple options available.
#' @param min.sz The genesets contain the minimum number of genes, where the default parameter is 1.
#' @param max.sz The genesets contain the maximum number of genes, where the default parameter is 10000.
#' @param p_combine_method The choice of method for integrating p-values.
#' @author Dingkang Zhao
#' @examples
#' @return
#' @export
DiffMetaScreener=function(expression_accession_vector,
                      comparisons,
                      geneSets_gmt,
                      enrichment_method=c("gsva_t","gsva_limma","gsva_anova","gsva_wilcoxon","gsva_permutation","gsva_kruskal",
                                         "ssgsea_t","ssgsea_limma","ssgsea_anova","ssgsea_wilcoxon","ssgsea_permutation","ssgsea_kruskal",
                                         "zscore_t","zscore_limma","zscore_anova","zscore_wilcoxon","zscore_permutation","zscore_kruskal",
                                         "plage_t","plage_limma","plage_anova","plage_wilcoxon","plage_permutation","plage_kruskal",
                                         "pca_t","pca_limma","pca_anova","pca_wilcoxon","pca_permutation","pca_kruskal",
                                         "aucell_t","aucell_limma","aucell_anova","aucell_wilcoxon","aucell_permutation","aucell_kruskal",
                                         "ucell_t","ucell_limma","ucell_anova","ucell_wilcoxon","ucell_permutation","ucell_kruskal",
                                         "singscore_t","singscore_limma","singscore_anova","singscore_wilcoxon","singscore_permutation","singscore_kruskal",
                                         "median_t","median_limma","median_anova","median_wilcoxon","median_permutation","median_kruskal",
                                         "t_fgsea","limma_fgsea","anova_fgsea","wilcoxon_fgsea","permutation_fgsea","kruskal_fgsea",
                                         "t_ora","limma_ora","anova_ora","wilcoxon_ora","permutation_ora","kruskal_ora",
                                         "consensus_t","consensus_limma","consensus_anova","consensus_wilcoxon","consensus_permutation","consensus_kruskal",
                                         "mdt_t","mdt_limma","mdt_anova","mdt_wilcoxon","mdt_permutation","mdt_kruskal",
                                         "mlm_t","mlm_limma","mlm_anova","mlm_wilcoxon","mlm_permutation","mlm_kruskal",
                                         "udt_t","udt_limma","udt_anova","udt_wilcoxon","udt_permutation","udt_kruskal",
                                         "ulm_t","ulm_limma","ulm_anova","ulm_wilcoxon","ulm_permutation","ulm_kruskal",
                                         "viper_t","viper_limma","viper_anova","viper_wilcoxon","viper_permutation","viper_kruskal",
                                         "wmean_t","wmean_limma","wmean_anova","wmean_wilcoxon","wmean_permutation","wmean_kruskal",
                                         "norm_wmean_t","norm_wmean_limma","norm_wmean_anova","norm_wmean_wilcoxon","norm_wmean_permutation","norm_wmean_kruskal",
                                         "corr_wmean_t","corr_wmean_limma","corr_wmean_anova","corr_wmean_wilcoxon","corr_wmean_permutation","corr_wmean_kruskal",
                                         "wsum_t","wsum_limma","wsum_anova","wsum_wilcoxon","wsum_permutation","wsum_kruskal",
                                         "norm_wsum_t","norm_wsum_limma","norm_wsum_anova","norm_wsum_wilcoxon","norm_wsum_permutation","norm_wsum_kruskal",
                                         "corr_wsum_t","corr_wsum_limma","corr_wsum_anova","corr_wsum_wilcoxon","corr_wsum_permutation","corr_wsum_kruskal"),
                       min.sz=1,
                       max.sz=10000,
                       p_combine_method=c("geometric_mean","invt",
                                          "invchisq","logitp","cct",
                                          "meanp","meanz","sumlog","sumz",
                                          "sump","votep","wilkinsonp")
){
  print("Start loading R packages.")
  library(singscore)
  library(UCell)
  library(decoupleR)
  library(AUCell)
  library(metap)
  library(IOBR)
  library(fgsea)
  library(plyr)
  library(reshape2)
  library(GSEABase)
  library(GSVA)
  library(doBy)
  library(clusterProfiler)
  library(limma)
  library(dplyr)
  library(tidyverse)
  library(devtools)
  library(stringr)
  library(GSA)
  library(coin)
  library(viper)
  ##########################################CTT: a p-value integration method
  #Liu Y, Xie J. Cauchy combination test: a powerful test with analytic p-value calculation under arbitrary dependency structures. J Am Stat Assoc. 2020;115(529):393-402. doi: 10.1080/01621459.2018.1554485. Epub 2019 Apr 25. PMID: 33012899; PMCID: PMC7531765.
  #Liu Y, Chen S, Li Z, Morrison AC, Boerwinkle E, Lin X. ACAT: A Fast and Powerful p Value Combination Method for Rare-Variant Analysis in Sequencing Studies. Am J Hum Genet. 2019 Mar 7;104(3):410-421. doi: 10.1016/j.ajhg.2019.01.002. PMID: 30849328; PMCID: PMC6407498.
  CCT = function(pvals, weights=NULL){
    if(sum(is.na(pvals)) > 0){
      stop("Cannot have NAs in the p-values!")
    }
    if((sum(pvals<0) + sum(pvals>1)) > 0){
      stop("All p-values must be between 0 and 1!")
    }
    is.zero = (sum(pvals==0)>=1)
    is.one = (sum(pvals==1)>=1)
    if(is.zero && is.one){
      zero_indices = which(pvals == 0)
      one_indices = which(pvals == 1)
      remove_count = min(length(zero_indices), length(one_indices))
      remove_zero = zero_indices[1:remove_count]
      remove_one = one_indices[1:remove_count]
      remove_indices = c(remove_zero, remove_one)
      pvals = pvals[-remove_indices]
      if (length(pvals) == 0){
        return(0)
      }
      if (!is.null(weights)) {
        weights = weights[-remove_indices]
      }
      warning("Removed ", remove_count, " p-value(s) of 0 and ", remove_count, " p-value(s) of 1 to avoid conflict.")
    }
    is.zero = (sum(pvals==0)>=1)
    is.one = (sum(pvals==1)>=1)
    if(is.zero){
      return(0)
    }
    if(is.one){
      warning("There are p-values that are exactly 1!")
      return(1)
    }
    if(is.null(weights)){
      weights = rep(1/length(pvals),length(pvals))
    }else if(length(weights)!=length(pvals)){
      stop("The length of weights should be the same as that of the p-values!")
    }else if(sum(weights < 0) > 0){
      stop("All the weights must be positive!")
    }else{
      weights = weights/sum(weights)
    }
    is.small = (pvals < 1e-16)
    if (sum(is.small) == 0){
      cct.stat = sum(weights*tan((0.5-pvals)*pi))
    }else{
      cct.stat = sum((weights[is.small]/pvals[is.small])/pi)
      cct.stat = cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
    }
    if(cct.stat > 1e+15){
      pval = (1/cct.stat)/pi
    }else{
      pval = 1-pcauchy(cct.stat)
    }
    return(pval)
  }
  ##########################################Median: The median of the genes contained in a gene set goes to measure the activity level of the gene set
  Median_evaluation=function(exp,geneSets){
    res=data.frame()
    for(num in 1:length(geneSets)){
      genes=geneSets[[num]]@geneIds
      genes_exp=exp[intersect(rownames(exp),genes),]
      genes_exp=as.data.frame(t(apply(genes_exp, 2, median)))
      rownames(genes_exp)=geneSets[[num]]@setName
      res=rbind(res,genes_exp)
    }
    res=na.omit(res)
    return(res)
  }
  ##########################################pvalue_integration_strategy: Parallel integration of multiple P-value integration methods.
  pvalue_integration_strategy=function(matrix,p_combine_method){
    if(length(colnames(matrix))>1){
      if(p_combine_method %in% c("invchisq")){
        p_combine_list=apply(matrix,1,function(row) invchisq(p=na.omit(row), k=length(na.omit(row))))
        p_combine=c()
        for(pnum in 1:length(p_combine_list)){
          p_combine_list_single=p_combine_list[[pnum]]$p
          names(p_combine_list_single)=names(p_combine_list)[pnum]
          p_combine=c(p_combine,p_combine_list_single)
        }
        p_combine=as.data.frame(p_combine)
      }else if(p_combine_method=="geometric_mean"){
        p_combine=apply(matrix,1,function(row) exp(mean(log(row),na.rm=T)))
        p_combine=as.data.frame(p_combine)
      }else if(p_combine_method=="invt"){
        p_combine_list=apply(matrix,1,function(row) invt(p=na.omit(row), k=length(na.omit(row))))
        p_combine=c()
        for(pnum in 1:length(p_combine_list)){
          p_combine_list_single=p_combine_list[[pnum]]$p
          names(p_combine_list_single)=names(p_combine_list)[pnum]
          p_combine=c(p_combine,p_combine_list_single)
        }
        p_combine=as.data.frame(p_combine)
      }else if(p_combine_method=="cct"){
        p_combine=apply(matrix,1,function(row) CCT(p=na.omit(row)))
        p_combine=as.data.frame(p_combine)
      }else if(p_combine_method=="logitp"){
        p_combine=apply(matrix,1,function(row) as.numeric(logitp(na.omit(row))$p))
        p_combine=as.data.frame(p_combine)
      }else if(p_combine_method=="meanp"){
        p_combine=apply(matrix,1,function(row) as.numeric(meanp(na.omit(row))$p))
        p_combine=as.data.frame(p_combine)
      }else if(p_combine_method=="meanz"){
        p_combine=apply(matrix,1,function(row) as.numeric(meanz(na.omit(row))$p))
        p_combine=as.data.frame(p_combine)
      }else if(p_combine_method=="sumlog"){
        p_combine=apply(matrix,1,function(row) as.numeric(sumlog(na.omit(row))$p))
        p_combine=as.data.frame(p_combine)
      }else if(p_combine_method=="sumz"){
        p_combine=apply(matrix,1,function(row) as.numeric(sumz(na.omit(row))$p))
        p_combine=as.data.frame(p_combine)
      }else if(p_combine_method=="sump"){
        p_combine=apply(matrix,1,function(row) as.numeric(sump(na.omit(row))$p))
        p_combine=as.data.frame(p_combine)
      }else if(p_combine_method=="votep"){
        p_combine=apply(matrix,1,function(row) as.numeric(votep(na.omit(row))$p))
        p_combine=as.data.frame(p_combine)
      }else if(p_combine_method=="wilkinsonp"){
        p_combine=apply(matrix,1,function(row) as.numeric(wilkinsonp(na.omit(row))$p))
        p_combine=as.data.frame(p_combine)
      }
    }else{
      p_combine=matrix
    }
    p_combine[is.na(p_combine)]=1
    return(p_combine)
  }
  ##########################################################################Load T-test function
  Get_T_test=function(geneset=geneset,group=group){
    result=NULL
    for (n in 1:nrow(geneset)) {
      geneset_n=data.frame(t(geneset[n,subset(group, group %in% c("activation","inhibition"))$Tag]))
      geneset_id=names(geneset_n)[1]
      names(geneset_n)[1]='geneset'
      geneset_n$Tag=rownames(geneset_n)
      geneset_n=merge(geneset_n, group, by = 'Tag', all.x = TRUE)
      geneset_n$group=factor(geneset_n$group,levels = c("activation","inhibition"))
      geneset_n_a=geneset_n[geneset_n$group=="activation",]
      geneset_n_i=geneset_n[geneset_n$group=="inhibition",]
      is_zero_variance=function(x) var(x) == 0
      if (is_zero_variance(geneset_n_a$geneset)&&is_zero_variance(geneset_n_i$geneset)) {
        next
      }
      P.Value=t.test(geneset~group, geneset_n,alternative="two.sided")$p.value
      if (!is.na(P.Value)) {
        stat=summaryBy(geneset~group, geneset_n, FUN = c(mean, median))
        result=rbind(result, c(geneset_id, as.character(stat[1,1]), stat[1,2], stat[1,3], as.character(stat[2,1]), stat[2,2], stat[2,3], P.Value))
      }
    }
    result=data.frame(result)
    names(result)=c('geneset_id', 'group1', 'mean1', 'median1', 'group2', 'mean2', 'median2', 'P.Value')
    result$mean1=as.numeric(result$mean1)
    result$mean2=as.numeric(result$mean2)
    result=result %>% mutate(logFC=mean1-mean2)
    result=as.data.frame(result)
    result$P.Value=as.numeric(result$P.Value)
    result$logFC=as.numeric(result$logFC)
    return(result) 
  }
  ##########################################################################Load wilcox-test function
  Get_wilcox_test=function(geneset=geneset,group=group){
    result=NULL
    for (n in 1:nrow(geneset)) {
      geneset_n=data.frame(t(geneset[n,subset(group, group %in% c("activation","inhibition"))$Tag]))
      geneset_id=names(geneset_n)[1]
      names(geneset_n)[1]='geneset'
      geneset_n$Tag=rownames(geneset_n)
      geneset_n=merge(geneset_n, group, by = 'Tag', all.x = TRUE)
      geneset_n$group=factor(geneset_n$group,levels = c("activation","inhibition"))
      geneset_n_a=geneset_n[geneset_n$group=="activation",]
      geneset_n_i=geneset_n[geneset_n$group=="inhibition",]
      is_zero_variance=function(x) var(x) == 0
      if (is_zero_variance(geneset_n_a$geneset)&&is_zero_variance(geneset_n_i$geneset)) {
        next
      }
      P.Value=wilcox.test(geneset~group, geneset_n,alternative="two.sided")$p.value
      if (!is.na(P.Value)) {
        stat=summaryBy(geneset~group, geneset_n, FUN = c(mean, median))
        result=rbind(result, c(geneset_id, as.character(stat[1,1]), stat[1,2], stat[1,3], as.character(stat[2,1]), stat[2,2], stat[2,3], P.Value))
      }
    }
    result=data.frame(result)
    names(result)=c('geneset_id', 'group1', 'mean1', 'median1', 'group2', 'mean2', 'median2', 'P.Value')
    result$mean1=as.numeric(result$mean1)
    result$mean2=as.numeric(result$mean2)
    result=result %>% mutate(logFC=mean1-mean2)
    result=as.data.frame(result)
    result$P.Value=as.numeric(result$P.Value)
    result$logFC=as.numeric(result$logFC)
    return(result)
  }
  ##########################################################################Load kruskal-test function
  Get_kruskal_test=function(geneset=geneset,group=group){
    result=NULL
    for (n in 1:nrow(geneset)) {
      geneset_n=data.frame(t(geneset[n,subset(group, group %in% c("activation","inhibition"))$Tag]))
      geneset_id=names(geneset_n)[1]
      names(geneset_n)[1]='geneset'
      geneset_n$Tag=rownames(geneset_n)
      geneset_n=merge(geneset_n, group, by = 'Tag', all.x = TRUE)
      geneset_n$group=factor(geneset_n$group,levels = c("activation","inhibition"))
      geneset_n_a=geneset_n[geneset_n$group=="activation",]
      geneset_n_i=geneset_n[geneset_n$group=="inhibition",]
      is_zero_variance=function(x) var(x) == 0
      if (is_zero_variance(geneset_n_a$geneset)&&is_zero_variance(geneset_n_i$geneset)) {
        next
      }
      P.Value=kruskal.test(geneset~group, geneset_n)$p.value
      if (!is.na(P.Value)) {
        stat=summaryBy(geneset~group, geneset_n, FUN = c(mean, median))
        result=rbind(result, c(geneset_id, as.character(stat[1,1]), stat[1,2], stat[1,3], as.character(stat[2,1]), stat[2,2], stat[2,3], P.Value))
      }
    }
    result=data.frame(result)
    names(result)=c('geneset_id', 'group1', 'mean1', 'median1', 'group2', 'mean2', 'median2', 'P.Value')
    result$mean1=as.numeric(result$mean1)
    result$mean2=as.numeric(result$mean2)
    result=result %>% mutate(logFC=mean1-mean2)
    result=as.data.frame(result)
    result$P.Value=as.numeric(result$P.Value)
    result$logFC=as.numeric(result$logFC)
    return(result)
  }
  ##########################################################################Load limma function
  Get_limma=function(geneset=geneset,group=group){
    design=model.matrix(~0+factor(group$group))
    colnames(design)=levels(factor(group$group))
    rownames(design)=colnames(geneset)
    compare=makeContrasts(activation-inhibition, levels=design)
    fit=lmFit(geneset, design)
    fit=contrasts.fit(fit, compare)
    fit=eBayes(fit)
    result=topTable(fit, coef=1, number=Inf, sort.by="none")
    result=na.omit(result)
    result$geneset_id=rownames(geneset)
    return(result)
  }
  ##########################################################################Load one-way ANOVA function
  Get_anova=function(geneset=geneset,group=group){
    result=NULL
    for (n in 1:nrow(geneset)) {
      geneset_n=data.frame(t(geneset[n,subset(group, group %in% c("activation","inhibition"))$Tag]))
      geneset_id=names(geneset_n)[1]
      names(geneset_n)[1]='geneset'
      geneset_n$Tag=rownames(geneset_n)
      geneset_n=merge(geneset_n, group, by = 'Tag', all.x = TRUE)
      geneset_n$group=factor(geneset_n$group,levels = c("activation","inhibition"))
      geneset_n_a=geneset_n[geneset_n$group=="activation",]
      geneset_n_i=geneset_n[geneset_n$group=="inhibition",]
      is_zero_variance=function(x) var(x) == 0
      if (is_zero_variance(geneset_n_a$geneset)&&is_zero_variance(geneset_n_i$geneset)) {
        next
      }
      P.Value=aov(geneset~group, geneset_n)
      P.Value=summary(P.Value)[[1]]
      P.Value=P.Value[1,5]
      if (!is.na(P.Value)) {
        stat=summaryBy(geneset~group, geneset_n, FUN = c(mean, median))
        result=rbind(result, c(geneset_id, as.character(stat[1,1]), stat[1,2], stat[1,3], as.character(stat[2,1]), stat[2,2], stat[2,3], P.Value))
      }
    }
    result=data.frame(result)
    names(result)=c('geneset_id', 'group1', 'mean1', 'median1', 'group2', 'mean2', 'median2', 'P.Value')
    result$mean1=as.numeric(result$mean1)
    result$mean2=as.numeric(result$mean2)
    result=result %>% mutate(logFC=mean1-mean2)
    result=as.data.frame(result)
    result$P.Value=as.numeric(result$P.Value)
    result$logFC=as.numeric(result$logFC)
    return(result)
  }
  ##########################################################################Load permutation-test function
  Get_permutation_test=function(geneset=geneset,group=group){
    geneset_group=merge(as.data.frame(t(geneset)),group,by.x="row.names",by.y="Tag")
    rownames(geneset_group)=geneset_group[,1]
    geneset_group=geneset_group[,-1]
    GS=colnames(geneset_group)
    GS=GS[GS!="group"]
    result=NULL
    for( n in 1:length(GS)){
      geneset_group_GS=geneset_group[,c(GS[n],"group")]
      geneset_group_GS$group=factor(geneset_group_GS$group,levels = c("activation","inhibition"))
      colnames(geneset_group_GS)[1]="geneset"
      geneset_n_a=geneset_group_GS[geneset_group_GS$group=="activation",]
      geneset_n_i=geneset_group_GS[geneset_group_GS$group=="inhibition",]
      is_zero_variance=function(x) var(x) == 0
      if (is_zero_variance(geneset_n_a$geneset)&&is_zero_variance(geneset_n_i$geneset)) {
        next
      }
      result_GS=oneway_test(geneset ~ group,data = geneset_group_GS,distribution = "asymptotic",alternative="two.sided")
      P.Value=pvalue(result_GS)
      mean1=geneset_group_GS[geneset_group_GS$group=="activation",]
      mean2=geneset_group_GS[geneset_group_GS$group=="inhibition",]
      mean1=mean(mean1[, 1], na.rm = TRUE)
      mean2=mean(mean2[, 1], na.rm = TRUE)
      result=rbind(result, c(GS[n],mean1,mean2,P.Value))
      
    }
    result=data.frame(result)
    names(result)=c('geneset_id',"mean1","mean2",'P.Value')
    result$mean1=as.numeric(result$mean1)
    result$mean2=as.numeric(result$mean2)
    result=result %>% mutate(logFC=mean1-mean2)
    result=as.data.frame(result)
    result$P.Value=as.numeric(result$P.Value)
    return(result)
  }
  ##########################################
  vector=expression_accession_vector
  print("Start using multiple enrichment analysis methods to cycle through the calculation of enrichment scores for each gene set in each dataset.")
  for (i in 1:length(vector)){
    print("Start loading expression data and grouping information.")
    print(paste0(vector[i]))
    exp=get(vector[i])
    dimnames=list(rownames(exp),colnames(exp))
    exp=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
    exp=as.data.frame(exp)
    constant_cols=sapply(exp, function(x) var(x) == 0)
    exp=exp[, !constant_cols]
    constant_rows=apply(exp, 1, function(x) length(unique(x)) == 1)
    exp=exp[!constant_rows, ]
    group=get(paste0(vector[i],"_G"))
    group=as.data.frame(group)
    Group_HL=comparisons
    rownames(Group_HL)=Group_HL$Accession
    Group_HL=Group_HL[,-1]
    group_HL_selected=Group_HL[paste0(vector[i]),,drop=F]
    group$group=ifelse(group$group==group_HL_selected$Group1,group_HL_selected$Group1_Status,group_HL_selected$Group2_Status)
    geneSets=geneSets_gmt
    label_genesets=names(geneSets)
    ##########################################################################Selection of the method of calculation
    ##########################################gsva
    methods_to_check=c("gsva_t", "gsva_limma", "gsva_anova", "gsva_wilcoxon", "gsva_permutation", "gsva_kruskal")
    if (any(methods_to_check %in% enrichment_method)){
      geneset=gsvaParam(exprData=as.matrix(exp),geneSets=geneSets,kcdf="Gaussian",minSize=min.sz,maxSize=max.sz)
      geneset=gsva(geneset)
      geneset=as.data.frame(geneset)
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("gsva_t" %in% enrichment_method){
        print("gsva_t will start.")
        result=Get_T_test(geneset=geneset,group=group)
        result_gsva_t=result[,c("geneset_id","P.Value","logFC")]
        result_gsva_t$two=TRUE
        result_gsva_t$invert_A=ifelse(result_gsva_t$logFC<=0,TRUE,FALSE)
        result_gsva_t$invert_I=ifelse(result_gsva_t$logFC>0,TRUE,FALSE)
        result_gsva_t$P.Value.A=two2one(result_gsva_t$P.Value, two = as.logical(result_gsva_t$two), invert = as.logical(result_gsva_t$invert_A))
        result_gsva_t$P.Value.I=two2one(result_gsva_t$P.Value, two = as.logical(result_gsva_t$two), invert = as.logical(result_gsva_t$invert_I))
        result_gsva_t_A=result_gsva_t[,c("geneset_id","P.Value.A")]
        result_gsva_t_I=result_gsva_t[,c("geneset_id","P.Value.I")]
        colnames(result_gsva_t_A)[2]=paste0(vector[i],".gsva_t")
        colnames(result_gsva_t_I)[2]=paste0(vector[i],".gsva_t")
        rm(result_gsva_t)
      }
      if("gsva_wilcoxon" %in% enrichment_method){
        print("gsva_wilcoxon will start.")
        result=Get_wilcox_test(geneset=geneset,group=group)
        result_gsva_wilcoxon=result[,c("geneset_id","P.Value","logFC")]
        result_gsva_wilcoxon$two=TRUE
        result_gsva_wilcoxon$invert_A=ifelse(result_gsva_wilcoxon$logFC<=0,TRUE,FALSE)
        result_gsva_wilcoxon$invert_I=ifelse(result_gsva_wilcoxon$logFC>0,TRUE,FALSE)
        result_gsva_wilcoxon$P.Value.A=two2one(result_gsva_wilcoxon$P.Value, two = as.logical(result_gsva_wilcoxon$two), invert = as.logical(result_gsva_wilcoxon$invert_A))
        result_gsva_wilcoxon$P.Value.I=two2one(result_gsva_wilcoxon$P.Value, two = as.logical(result_gsva_wilcoxon$two), invert = as.logical(result_gsva_wilcoxon$invert_I))
        result_gsva_wilcoxon_A=result_gsva_wilcoxon[,c("geneset_id","P.Value.A")]
        result_gsva_wilcoxon_I=result_gsva_wilcoxon[,c("geneset_id","P.Value.I")]
        colnames(result_gsva_wilcoxon_A)[2]=paste0(vector[i],".gsva_wilcoxon")
        colnames(result_gsva_wilcoxon_I)[2]=paste0(vector[i],".gsva_wilcoxon")
        rm(result_gsva_wilcoxon)
      }
      if("gsva_kruskal" %in% enrichment_method){
        print("gsva_kruskal will start.")
        result=Get_kruskal_test(geneset=geneset,group=group)
        result_gsva_kruskal=result[,c("geneset_id","P.Value","logFC")]
        result_gsva_kruskal$two=TRUE
        result_gsva_kruskal$invert_A=ifelse(result_gsva_kruskal$logFC<=0,TRUE,FALSE)
        result_gsva_kruskal$invert_I=ifelse(result_gsva_kruskal$logFC>0,TRUE,FALSE)
        result_gsva_kruskal$P.Value.A=two2one(result_gsva_kruskal$P.Value, two = as.logical(result_gsva_kruskal$two), invert = as.logical(result_gsva_kruskal$invert_A))
        result_gsva_kruskal$P.Value.I=two2one(result_gsva_kruskal$P.Value, two = as.logical(result_gsva_kruskal$two), invert = as.logical(result_gsva_kruskal$invert_I))
        result_gsva_kruskal_A=result_gsva_kruskal[,c("geneset_id","P.Value.A")]
        result_gsva_kruskal_I=result_gsva_kruskal[,c("geneset_id","P.Value.I")]
        colnames(result_gsva_kruskal_A)[2]=paste0(vector[i],".gsva_kruskal")
        colnames(result_gsva_kruskal_I)[2]=paste0(vector[i],".gsva_kruskal")
        rm(result_gsva_kruskal)
      }
      if("gsva_permutation" %in% enrichment_method){
        print("gsva_permutation will start.")
        result=Get_permutation_test(geneset=geneset,group=group)
        result_gsva_permutation=result[,c("geneset_id","P.Value","logFC")]
        result_gsva_permutation$two=TRUE
        result_gsva_permutation$invert_A=ifelse(result_gsva_permutation$logFC<=0,TRUE,FALSE)
        result_gsva_permutation$invert_I=ifelse(result_gsva_permutation$logFC>0,TRUE,FALSE)
        result_gsva_permutation=na.omit(result_gsva_permutation)
        result_gsva_permutation$P.Value.A=two2one(result_gsva_permutation$P.Value, two = as.logical(result_gsva_permutation$two), invert = as.logical(result_gsva_permutation$invert_A))
        result_gsva_permutation$P.Value.I=two2one(result_gsva_permutation$P.Value, two = as.logical(result_gsva_permutation$two), invert = as.logical(result_gsva_permutation$invert_I))
        result_gsva_permutation_A=result_gsva_permutation[,c("geneset_id","P.Value.A")]
        result_gsva_permutation_I=result_gsva_permutation[,c("geneset_id","P.Value.I")]
        colnames(result_gsva_permutation_A)[2]=paste0(vector[i],".gsva_permutation")
        colnames(result_gsva_permutation_I)[2]=paste0(vector[i],".gsva_permutation")
        rm(result_gsva_permutation)
      }
      if("gsva_limma" %in% enrichment_method){
        print("gsva_limma will start.")
        result=Get_limma(geneset=geneset,group=group)
        result_gsva_limma=result[,c("geneset_id","P.Value","logFC")]
        result_gsva_limma$two=TRUE
        result_gsva_limma$invert_A=ifelse(result_gsva_limma$logFC<=0,TRUE,FALSE)
        result_gsva_limma$invert_I=ifelse(result_gsva_limma$logFC>0,TRUE,FALSE)
        result_gsva_limma$P.Value.A=two2one(result_gsva_limma$P.Value, two = as.logical(result_gsva_limma$two), invert = as.logical(result_gsva_limma$invert_A))
        result_gsva_limma$P.Value.I=two2one(result_gsva_limma$P.Value, two = as.logical(result_gsva_limma$two), invert = as.logical(result_gsva_limma$invert_I))
        result_gsva_limma_A=result_gsva_limma[,c("geneset_id","P.Value.A")]
        result_gsva_limma_I=result_gsva_limma[,c("geneset_id","P.Value.I")]
        colnames(result_gsva_limma_A)[2]=paste0(vector[i],".gsva_limma")
        colnames(result_gsva_limma_I)[2]=paste0(vector[i],".gsva_limma")
        rm(result_gsva_limma)
      }
      if("gsva_anova" %in% enrichment_method){
        print("gsva_anova will start.")
        result=Get_anova(geneset=geneset,group=group)
        result_gsva_anova=result[,c("geneset_id","P.Value","logFC")]
        result_gsva_anova$two=TRUE
        result_gsva_anova$invert_A=ifelse(result_gsva_anova$logFC<=0,TRUE,FALSE)
        result_gsva_anova$invert_I=ifelse(result_gsva_anova$logFC>0,TRUE,FALSE)
        result_gsva_anova$P.Value.A=two2one(result_gsva_anova$P.Value, two = as.logical(result_gsva_anova$two), invert = as.logical(result_gsva_anova$invert_A))
        result_gsva_anova$P.Value.I=two2one(result_gsva_anova$P.Value, two = as.logical(result_gsva_anova$two), invert = as.logical(result_gsva_anova$invert_I))
        result_gsva_anova_A=result_gsva_anova[,c("geneset_id","P.Value.A")]
        result_gsva_anova_I=result_gsva_anova[,c("geneset_id","P.Value.I")]
        colnames(result_gsva_anova_A)[2]=paste0(vector[i],".gsva_anova")
        colnames(result_gsva_anova_I)[2]=paste0(vector[i],".gsva_anova")
        rm(result_gsva_anova)
      }
    } else {
      print("Not using gsva.")
    }
    ##########################################ssgsea
    methods_to_check=c("ssgsea_t", "ssgsea_limma", "ssgsea_anova", "ssgsea_wilcoxon", "ssgsea_permutation", "ssgsea_kruskal")
    if (any(methods_to_check %in% enrichment_method)){
      geneset=ssgseaParam(exprData=as.matrix(exp),geneSets=geneSets,minSize=min.sz,maxSize=max.sz)
      geneset=gsva(geneset)
      geneset=as.data.frame(geneset)
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("ssgsea_t" %in% enrichment_method){
        print("ssgsea_t will start.")
        result=Get_T_test(geneset=geneset,group=group)
        result_ssgsea_t=result[,c("geneset_id","P.Value","logFC")]
        result_ssgsea_t$two=TRUE
        result_ssgsea_t$invert_A=ifelse(result_ssgsea_t$logFC<=0,TRUE,FALSE)
        result_ssgsea_t$invert_I=ifelse(result_ssgsea_t$logFC>0,TRUE,FALSE)
        result_ssgsea_t$P.Value.A=two2one(result_ssgsea_t$P.Value, two = as.logical(result_ssgsea_t$two), invert = as.logical(result_ssgsea_t$invert_A))
        result_ssgsea_t$P.Value.I=two2one(result_ssgsea_t$P.Value, two = as.logical(result_ssgsea_t$two), invert = as.logical(result_ssgsea_t$invert_I))
        result_ssgsea_t_A=result_ssgsea_t[,c("geneset_id","P.Value.A")]
        result_ssgsea_t_I=result_ssgsea_t[,c("geneset_id","P.Value.I")]
        colnames(result_ssgsea_t_A)[2]=paste0(vector[i],".ssgsea_t")
        colnames(result_ssgsea_t_I)[2]=paste0(vector[i],".ssgsea_t")
        rm(result_ssgsea_t)
      }
      if("ssgsea_wilcoxon" %in% enrichment_method){
        print("ssgsea_wilcoxon will start.")
        result=Get_wilcox_test(geneset=geneset,group=group)
        result_ssgsea_wilcoxon=result[,c("geneset_id","P.Value","logFC")]
        result_ssgsea_wilcoxon$two=TRUE
        result_ssgsea_wilcoxon$invert_A=ifelse(result_ssgsea_wilcoxon$logFC<=0,TRUE,FALSE)
        result_ssgsea_wilcoxon$invert_I=ifelse(result_ssgsea_wilcoxon$logFC>0,TRUE,FALSE)
        result_ssgsea_wilcoxon$P.Value.A=two2one(result_ssgsea_wilcoxon$P.Value, two = as.logical(result_ssgsea_wilcoxon$two), invert = as.logical(result_ssgsea_wilcoxon$invert_A))
        result_ssgsea_wilcoxon$P.Value.I=two2one(result_ssgsea_wilcoxon$P.Value, two = as.logical(result_ssgsea_wilcoxon$two), invert = as.logical(result_ssgsea_wilcoxon$invert_I))
        result_ssgsea_wilcoxon_A=result_ssgsea_wilcoxon[,c("geneset_id","P.Value.A")]
        result_ssgsea_wilcoxon_I=result_ssgsea_wilcoxon[,c("geneset_id","P.Value.I")]
        colnames(result_ssgsea_wilcoxon_A)[2]=paste0(vector[i],".ssgsea_wilcoxon")
        colnames(result_ssgsea_wilcoxon_I)[2]=paste0(vector[i],".ssgsea_wilcoxon")
        rm(result_ssgsea_wilcoxon)
      }
      if("ssgsea_kruskal" %in% enrichment_method){
        print("ssgsea_kruskal will start.")
        result=Get_kruskal_test(geneset=geneset,group=group)
        result_ssgsea_kruskal=result[,c("geneset_id","P.Value","logFC")]
        result_ssgsea_kruskal$two=TRUE
        result_ssgsea_kruskal$invert_A=ifelse(result_ssgsea_kruskal$logFC<=0,TRUE,FALSE)
        result_ssgsea_kruskal$invert_I=ifelse(result_ssgsea_kruskal$logFC>0,TRUE,FALSE)
        result_ssgsea_kruskal$P.Value.A=two2one(result_ssgsea_kruskal$P.Value, two = as.logical(result_ssgsea_kruskal$two), invert = as.logical(result_ssgsea_kruskal$invert_A))
        result_ssgsea_kruskal$P.Value.I=two2one(result_ssgsea_kruskal$P.Value, two = as.logical(result_ssgsea_kruskal$two), invert = as.logical(result_ssgsea_kruskal$invert_I))
        result_ssgsea_kruskal_A=result_ssgsea_kruskal[,c("geneset_id","P.Value.A")]
        result_ssgsea_kruskal_I=result_ssgsea_kruskal[,c("geneset_id","P.Value.I")]
        colnames(result_ssgsea_kruskal_A)[2]=paste0(vector[i],".ssgsea_kruskal")
        colnames(result_ssgsea_kruskal_I)[2]=paste0(vector[i],".ssgsea_kruskal")
        rm(result_ssgsea_kruskal)
      }
      if("ssgsea_permutation" %in% enrichment_method){
        print("ssgsea_permutation will start.")
        result=Get_permutation_test(geneset=geneset,group=group)
        result_ssgsea_permutation=result[,c("geneset_id","P.Value","logFC")]
        result_ssgsea_permutation$two=TRUE
        result_ssgsea_permutation$invert_A=ifelse(result_ssgsea_permutation$logFC<=0,TRUE,FALSE)
        result_ssgsea_permutation$invert_I=ifelse(result_ssgsea_permutation$logFC>0,TRUE,FALSE)
        result_ssgsea_permutation=na.omit(result_ssgsea_permutation)
        result_ssgsea_permutation$P.Value.A=two2one(result_ssgsea_permutation$P.Value, two = as.logical(result_ssgsea_permutation$two), invert = as.logical(result_ssgsea_permutation$invert_A))
        result_ssgsea_permutation$P.Value.I=two2one(result_ssgsea_permutation$P.Value, two = as.logical(result_ssgsea_permutation$two), invert = as.logical(result_ssgsea_permutation$invert_I))
        result_ssgsea_permutation_A=result_ssgsea_permutation[,c("geneset_id","P.Value.A")]
        result_ssgsea_permutation_I=result_ssgsea_permutation[,c("geneset_id","P.Value.I")]
        colnames(result_ssgsea_permutation_A)[2]=paste0(vector[i],".ssgsea_permutation")
        colnames(result_ssgsea_permutation_I)[2]=paste0(vector[i],".ssgsea_permutation")
        rm(result_ssgsea_permutation)
      }
      if("ssgsea_limma" %in% enrichment_method){
        print("ssgsea_limma will start.")
        result=Get_limma(geneset=geneset,group=group)
        result_ssgsea_limma=result[,c("geneset_id","P.Value","logFC")]
        result_ssgsea_limma$two=TRUE
        result_ssgsea_limma$invert_A=ifelse(result_ssgsea_limma$logFC<=0,TRUE,FALSE)
        result_ssgsea_limma$invert_I=ifelse(result_ssgsea_limma$logFC>0,TRUE,FALSE)
        result_ssgsea_limma$P.Value.A=two2one(result_ssgsea_limma$P.Value, two = as.logical(result_ssgsea_limma$two), invert = as.logical(result_ssgsea_limma$invert_A))
        result_ssgsea_limma$P.Value.I=two2one(result_ssgsea_limma$P.Value, two = as.logical(result_ssgsea_limma$two), invert = as.logical(result_ssgsea_limma$invert_I))
        result_ssgsea_limma_A=result_ssgsea_limma[,c("geneset_id","P.Value.A")]
        result_ssgsea_limma_I=result_ssgsea_limma[,c("geneset_id","P.Value.I")]
        colnames(result_ssgsea_limma_A)[2]=paste0(vector[i],".ssgsea_limma")
        colnames(result_ssgsea_limma_I)[2]=paste0(vector[i],".ssgsea_limma")
        rm(result_ssgsea_limma)
      }
      if("ssgsea_anova" %in% enrichment_method){
        print("ssgsea_anova will start.")
        result=Get_anova(geneset=geneset,group=group)
        result_ssgsea_anova=result[,c("geneset_id","P.Value","logFC")]
        result_ssgsea_anova$two=TRUE
        result_ssgsea_anova$invert_A=ifelse(result_ssgsea_anova$logFC<=0,TRUE,FALSE)
        result_ssgsea_anova$invert_I=ifelse(result_ssgsea_anova$logFC>0,TRUE,FALSE)
        result_ssgsea_anova$P.Value.A=two2one(result_ssgsea_anova$P.Value, two = as.logical(result_ssgsea_anova$two), invert = as.logical(result_ssgsea_anova$invert_A))
        result_ssgsea_anova$P.Value.I=two2one(result_ssgsea_anova$P.Value, two = as.logical(result_ssgsea_anova$two), invert = as.logical(result_ssgsea_anova$invert_I))
        result_ssgsea_anova_A=result_ssgsea_anova[,c("geneset_id","P.Value.A")]
        result_ssgsea_anova_I=result_ssgsea_anova[,c("geneset_id","P.Value.I")]
        colnames(result_ssgsea_anova_A)[2]=paste0(vector[i],".ssgsea_anova")
        colnames(result_ssgsea_anova_I)[2]=paste0(vector[i],".ssgsea_anova")
        rm(result_ssgsea_anova)
      }
    } else {
      print("Not using ssgsea.")
    }
    ##########################################zscore
    methods_to_check=c("zscore_t", "zscore_limma", "zscore_anova", "zscore_wilcoxon", "zscore_permutation", "zscore_kruskal")
    if (any(methods_to_check %in% enrichment_method)){
      geneset=zscoreParam(exprData=as.matrix(exp),geneSets=geneSets,minSize=min.sz,maxSize=max.sz)
      geneset=gsva(geneset)
      geneset=as.data.frame(geneset)
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("zscore_t" %in% enrichment_method){
        print("zscore_t will start.")
        result=Get_T_test(geneset=geneset,group=group)
        result_zscore_t=result[,c("geneset_id","P.Value","logFC")]
        result_zscore_t$two=TRUE
        result_zscore_t$invert_A=ifelse(result_zscore_t$logFC<=0,TRUE,FALSE)
        result_zscore_t$invert_I=ifelse(result_zscore_t$logFC>0,TRUE,FALSE)
        result_zscore_t$P.Value.A=two2one(result_zscore_t$P.Value, two = as.logical(result_zscore_t$two), invert = as.logical(result_zscore_t$invert_A))
        result_zscore_t$P.Value.I=two2one(result_zscore_t$P.Value, two = as.logical(result_zscore_t$two), invert = as.logical(result_zscore_t$invert_I))
        result_zscore_t_A=result_zscore_t[,c("geneset_id","P.Value.A")]
        result_zscore_t_I=result_zscore_t[,c("geneset_id","P.Value.I")]
        colnames(result_zscore_t_A)[2]=paste0(vector[i],".zscore_t")
        colnames(result_zscore_t_I)[2]=paste0(vector[i],".zscore_t")
        rm(result_zscore_t)
      }
      if("zscore_wilcoxon" %in% enrichment_method){
        print("zscore_wilcoxon will start.")
        result=Get_wilcox_test(geneset=geneset,group=group)
        result_zscore_wilcoxon=result[,c("geneset_id","P.Value","logFC")]
        result_zscore_wilcoxon$two=TRUE
        result_zscore_wilcoxon$invert_A=ifelse(result_zscore_wilcoxon$logFC<=0,TRUE,FALSE)
        result_zscore_wilcoxon$invert_I=ifelse(result_zscore_wilcoxon$logFC>0,TRUE,FALSE)
        result_zscore_wilcoxon$P.Value.A=two2one(result_zscore_wilcoxon$P.Value, two = as.logical(result_zscore_wilcoxon$two), invert = as.logical(result_zscore_wilcoxon$invert_A))
        result_zscore_wilcoxon$P.Value.I=two2one(result_zscore_wilcoxon$P.Value, two = as.logical(result_zscore_wilcoxon$two), invert = as.logical(result_zscore_wilcoxon$invert_I))
        result_zscore_wilcoxon_A=result_zscore_wilcoxon[,c("geneset_id","P.Value.A")]
        result_zscore_wilcoxon_I=result_zscore_wilcoxon[,c("geneset_id","P.Value.I")]
        colnames(result_zscore_wilcoxon_A)[2]=paste0(vector[i],".zscore_wilcoxon")
        colnames(result_zscore_wilcoxon_I)[2]=paste0(vector[i],".zscore_wilcoxon")
        rm(result_zscore_wilcoxon)
      }
      if("zscore_kruskal" %in% enrichment_method){
        print("zscore_kruskal will start.")
        result=Get_kruskal_test(geneset=geneset,group=group)
        result_zscore_kruskal=result[,c("geneset_id","P.Value","logFC")]
        result_zscore_kruskal$two=TRUE
        result_zscore_kruskal$invert_A=ifelse(result_zscore_kruskal$logFC<=0,TRUE,FALSE)
        result_zscore_kruskal$invert_I=ifelse(result_zscore_kruskal$logFC>0,TRUE,FALSE)
        result_zscore_kruskal$P.Value.A=two2one(result_zscore_kruskal$P.Value, two = as.logical(result_zscore_kruskal$two), invert = as.logical(result_zscore_kruskal$invert_A))
        result_zscore_kruskal$P.Value.I=two2one(result_zscore_kruskal$P.Value, two = as.logical(result_zscore_kruskal$two), invert = as.logical(result_zscore_kruskal$invert_I))
        result_zscore_kruskal_A=result_zscore_kruskal[,c("geneset_id","P.Value.A")]
        result_zscore_kruskal_I=result_zscore_kruskal[,c("geneset_id","P.Value.I")]
        colnames(result_zscore_kruskal_A)[2]=paste0(vector[i],".zscore_kruskal")
        colnames(result_zscore_kruskal_I)[2]=paste0(vector[i],".zscore_kruskal")
        rm(result_zscore_kruskal)
      }
      if("zscore_permutation" %in% enrichment_method){
        print("zscore_permutation will start.")
        result=Get_permutation_test(geneset=geneset,group=group)
        result_zscore_permutation=result[,c("geneset_id","P.Value","logFC")]
        result_zscore_permutation$two=TRUE
        result_zscore_permutation$invert_A=ifelse(result_zscore_permutation$logFC<=0,TRUE,FALSE)
        result_zscore_permutation$invert_I=ifelse(result_zscore_permutation$logFC>0,TRUE,FALSE)
        result_zscore_permutation=na.omit(result_zscore_permutation)
        result_zscore_permutation$P.Value.A=two2one(result_zscore_permutation$P.Value, two = as.logical(result_zscore_permutation$two), invert = as.logical(result_zscore_permutation$invert_A))
        result_zscore_permutation$P.Value.I=two2one(result_zscore_permutation$P.Value, two = as.logical(result_zscore_permutation$two), invert = as.logical(result_zscore_permutation$invert_I))
        result_zscore_permutation_A=result_zscore_permutation[,c("geneset_id","P.Value.A")]
        result_zscore_permutation_I=result_zscore_permutation[,c("geneset_id","P.Value.I")]
        colnames(result_zscore_permutation_A)[2]=paste0(vector[i],".zscore_permutation")
        colnames(result_zscore_permutation_I)[2]=paste0(vector[i],".zscore_permutation")
        rm(result_zscore_permutation)
      }
      if("zscore_limma" %in% enrichment_method){
        print("zscore_limma will start.")
        result=Get_limma(geneset=geneset,group=group)
        result_zscore_limma=result[,c("geneset_id","P.Value","logFC")]
        result_zscore_limma$two=TRUE
        result_zscore_limma$invert_A=ifelse(result_zscore_limma$logFC<=0,TRUE,FALSE)
        result_zscore_limma$invert_I=ifelse(result_zscore_limma$logFC>0,TRUE,FALSE)
        result_zscore_limma$P.Value.A=two2one(result_zscore_limma$P.Value, two = as.logical(result_zscore_limma$two), invert = as.logical(result_zscore_limma$invert_A))
        result_zscore_limma$P.Value.I=two2one(result_zscore_limma$P.Value, two = as.logical(result_zscore_limma$two), invert = as.logical(result_zscore_limma$invert_I))
        result_zscore_limma_A=result_zscore_limma[,c("geneset_id","P.Value.A")]
        result_zscore_limma_I=result_zscore_limma[,c("geneset_id","P.Value.I")]
        colnames(result_zscore_limma_A)[2]=paste0(vector[i],".zscore_limma")
        colnames(result_zscore_limma_I)[2]=paste0(vector[i],".zscore_limma")
        rm(result_zscore_limma)
      }
      if("zscore_anova" %in% enrichment_method){
        print("zscore_anova will start.")
        result=Get_anova(geneset=geneset,group=group)
        result_zscore_anova=result[,c("geneset_id","P.Value","logFC")]
        result_zscore_anova$two=TRUE
        result_zscore_anova$invert_A=ifelse(result_zscore_anova$logFC<=0,TRUE,FALSE)
        result_zscore_anova$invert_I=ifelse(result_zscore_anova$logFC>0,TRUE,FALSE)
        result_zscore_anova$P.Value.A=two2one(result_zscore_anova$P.Value, two = as.logical(result_zscore_anova$two), invert = as.logical(result_zscore_anova$invert_A))
        result_zscore_anova$P.Value.I=two2one(result_zscore_anova$P.Value, two = as.logical(result_zscore_anova$two), invert = as.logical(result_zscore_anova$invert_I))
        result_zscore_anova_A=result_zscore_anova[,c("geneset_id","P.Value.A")]
        result_zscore_anova_I=result_zscore_anova[,c("geneset_id","P.Value.I")]
        colnames(result_zscore_anova_A)[2]=paste0(vector[i],".zscore_anova")
        colnames(result_zscore_anova_I)[2]=paste0(vector[i],".zscore_anova")
        rm(result_zscore_anova)
      }
    } else {
      print("Not using zscore.")
    }
    ##########################################plage
    methods_to_check=c("plage_t", "plage_limma", "plage_anova", "plage_wilcoxon", "plage_permutation", "plage_kruskal")
    if (any(methods_to_check %in% enrichment_method)){
      geneset=plageParam(exprData=as.matrix(exp),geneSets=geneSets,minSize=min.sz,maxSize=max.sz)
      geneset=gsva(geneset)
      geneset=as.data.frame(geneset)
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("plage_t" %in% enrichment_method){
        print("plage_t will start.")
        result=Get_T_test(geneset=geneset,group=group)
        result_plage_t=result[,c("geneset_id","P.Value","logFC")]
        result_plage_t$two=TRUE
        result_plage_t$invert_A=ifelse(result_plage_t$logFC<=0,TRUE,FALSE)
        result_plage_t$invert_I=ifelse(result_plage_t$logFC>0,TRUE,FALSE)
        result_plage_t$P.Value.A=two2one(result_plage_t$P.Value, two = as.logical(result_plage_t$two), invert = as.logical(result_plage_t$invert_A))
        result_plage_t$P.Value.I=two2one(result_plage_t$P.Value, two = as.logical(result_plage_t$two), invert = as.logical(result_plage_t$invert_I))
        result_plage_t_A=result_plage_t[,c("geneset_id","P.Value.A")]
        result_plage_t_I=result_plage_t[,c("geneset_id","P.Value.I")]
        colnames(result_plage_t_A)[2]=paste0(vector[i],".plage_t")
        colnames(result_plage_t_I)[2]=paste0(vector[i],".plage_t")
        rm(result_plage_t)
      }
      if("plage_wilcoxon" %in% enrichment_method){
        print("plage_wilcoxon will start.")
        result=Get_wilcox_test(geneset=geneset,group=group)
        result_plage_wilcoxon=result[,c("geneset_id","P.Value","logFC")]
        result_plage_wilcoxon$two=TRUE
        result_plage_wilcoxon$invert_A=ifelse(result_plage_wilcoxon$logFC<=0,TRUE,FALSE)
        result_plage_wilcoxon$invert_I=ifelse(result_plage_wilcoxon$logFC>0,TRUE,FALSE)
        result_plage_wilcoxon$P.Value.A=two2one(result_plage_wilcoxon$P.Value, two = as.logical(result_plage_wilcoxon$two), invert = as.logical(result_plage_wilcoxon$invert_A))
        result_plage_wilcoxon$P.Value.I=two2one(result_plage_wilcoxon$P.Value, two = as.logical(result_plage_wilcoxon$two), invert = as.logical(result_plage_wilcoxon$invert_I))
        result_plage_wilcoxon_A=result_plage_wilcoxon[,c("geneset_id","P.Value.A")]
        result_plage_wilcoxon_I=result_plage_wilcoxon[,c("geneset_id","P.Value.I")]
        colnames(result_plage_wilcoxon_A)[2]=paste0(vector[i],".plage_wilcoxon")
        colnames(result_plage_wilcoxon_I)[2]=paste0(vector[i],".plage_wilcoxon")
        rm(result_plage_wilcoxon)
      }
      if("plage_kruskal" %in% enrichment_method){
        print("plage_kruskal will start.")
        result=Get_kruskal_test(geneset=geneset,group=group)
        result_plage_kruskal=result[,c("geneset_id","P.Value","logFC")]
        result_plage_kruskal$two=TRUE
        result_plage_kruskal$invert_A=ifelse(result_plage_kruskal$logFC<=0,TRUE,FALSE)
        result_plage_kruskal$invert_I=ifelse(result_plage_kruskal$logFC>0,TRUE,FALSE)
        result_plage_kruskal$P.Value.A=two2one(result_plage_kruskal$P.Value, two = as.logical(result_plage_kruskal$two), invert = as.logical(result_plage_kruskal$invert_A))
        result_plage_kruskal$P.Value.I=two2one(result_plage_kruskal$P.Value, two = as.logical(result_plage_kruskal$two), invert = as.logical(result_plage_kruskal$invert_I))
        result_plage_kruskal_A=result_plage_kruskal[,c("geneset_id","P.Value.A")]
        result_plage_kruskal_I=result_plage_kruskal[,c("geneset_id","P.Value.I")]
        colnames(result_plage_kruskal_A)[2]=paste0(vector[i],".plage_kruskal")
        colnames(result_plage_kruskal_I)[2]=paste0(vector[i],".plage_kruskal")
        rm(result_plage_kruskal)
      }
      if("plage_permutation" %in% enrichment_method){
        print("plage_permutation will start.")
        result=Get_permutation_test(geneset=geneset,group=group)
        result_plage_permutation=result[,c("geneset_id","P.Value","logFC")]
        result_plage_permutation$two=TRUE
        result_plage_permutation$invert_A=ifelse(result_plage_permutation$logFC<=0,TRUE,FALSE)
        result_plage_permutation$invert_I=ifelse(result_plage_permutation$logFC>0,TRUE,FALSE)
        result_plage_permutation=na.omit(result_plage_permutation)
        result_plage_permutation$P.Value.A=two2one(result_plage_permutation$P.Value, two = as.logical(result_plage_permutation$two), invert = as.logical(result_plage_permutation$invert_A))
        result_plage_permutation$P.Value.I=two2one(result_plage_permutation$P.Value, two = as.logical(result_plage_permutation$two), invert = as.logical(result_plage_permutation$invert_I))
        result_plage_permutation_A=result_plage_permutation[,c("geneset_id","P.Value.A")]
        result_plage_permutation_I=result_plage_permutation[,c("geneset_id","P.Value.I")]
        colnames(result_plage_permutation_A)[2]=paste0(vector[i],".plage_permutation")
        colnames(result_plage_permutation_I)[2]=paste0(vector[i],".plage_permutation")
        rm(result_plage_permutation)
      }
      if("plage_limma" %in% enrichment_method){
        print("plage_limma will start.")
        result=Get_limma(geneset=geneset,group=group)
        result_plage_limma=result[,c("geneset_id","P.Value","logFC")]
        result_plage_limma$two=TRUE
        result_plage_limma$invert_A=ifelse(result_plage_limma$logFC<=0,TRUE,FALSE)
        result_plage_limma$invert_I=ifelse(result_plage_limma$logFC>0,TRUE,FALSE)
        result_plage_limma$P.Value.A=two2one(result_plage_limma$P.Value, two = as.logical(result_plage_limma$two), invert = as.logical(result_plage_limma$invert_A))
        result_plage_limma$P.Value.I=two2one(result_plage_limma$P.Value, two = as.logical(result_plage_limma$two), invert = as.logical(result_plage_limma$invert_I))
        result_plage_limma_A=result_plage_limma[,c("geneset_id","P.Value.A")]
        result_plage_limma_I=result_plage_limma[,c("geneset_id","P.Value.I")]
        colnames(result_plage_limma_A)[2]=paste0(vector[i],".plage_limma")
        colnames(result_plage_limma_I)[2]=paste0(vector[i],".plage_limma")
        rm(result_plage_limma)
      }
      if("plage_anova" %in% enrichment_method){
        print("plage_anova will start.")
        result=Get_anova(geneset=geneset,group=group)
        result_plage_anova=result[,c("geneset_id","P.Value","logFC")]
        result_plage_anova$two=TRUE
        result_plage_anova$invert_A=ifelse(result_plage_anova$logFC<=0,TRUE,FALSE)
        result_plage_anova$invert_I=ifelse(result_plage_anova$logFC>0,TRUE,FALSE)
        result_plage_anova$P.Value.A=two2one(result_plage_anova$P.Value, two = as.logical(result_plage_anova$two), invert = as.logical(result_plage_anova$invert_A))
        result_plage_anova$P.Value.I=two2one(result_plage_anova$P.Value, two = as.logical(result_plage_anova$two), invert = as.logical(result_plage_anova$invert_I))
        result_plage_anova_A=result_plage_anova[,c("geneset_id","P.Value.A")]
        result_plage_anova_I=result_plage_anova[,c("geneset_id","P.Value.I")]
        colnames(result_plage_anova_A)[2]=paste0(vector[i],".plage_anova")
        colnames(result_plage_anova_I)[2]=paste0(vector[i],".plage_anova")
        rm(result_plage_anova)
      }
    } else {
      print("Not using plage.")
    }
    ##########################################pca
    methods_to_check=c("pca_t", "pca_limma", "pca_anova", "pca_wilcoxon", "pca_permutation", "pca_kruskal")
    if (any(methods_to_check %in% enrichment_method)){
      geneSets_list=lapply(geneSets, function(x){ x@geneIds })
      names(geneSets_list)=names(geneSets)
      geneset=calculate_sig_score(eset=as.matrix(exp),
                                  signature=geneSets_list,
                                  method="pca",
                                  mini_gene_count=min.sz)
      geneset=as.data.frame(geneset)
      rownames(geneset)=geneset[,1]
      geneset=geneset[,-1]
      geneset=as.data.frame(t(geneset))
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("pca_t" %in% enrichment_method){
        print("pca_t will start.")
        result=Get_T_test(geneset=geneset,group=group)
        result_pca_t=result[,c("geneset_id","P.Value","logFC")]
        result_pca_t$two=TRUE
        result_pca_t$invert_A=ifelse(result_pca_t$logFC<=0,TRUE,FALSE)
        result_pca_t$invert_I=ifelse(result_pca_t$logFC>0,TRUE,FALSE)
        result_pca_t$P.Value.A=two2one(result_pca_t$P.Value, two = as.logical(result_pca_t$two), invert = as.logical(result_pca_t$invert_A))
        result_pca_t$P.Value.I=two2one(result_pca_t$P.Value, two = as.logical(result_pca_t$two), invert = as.logical(result_pca_t$invert_I))
        result_pca_t_A=result_pca_t[,c("geneset_id","P.Value.A")]
        result_pca_t_I=result_pca_t[,c("geneset_id","P.Value.I")]
        colnames(result_pca_t_A)[2]=paste0(vector[i],".pca_t")
        colnames(result_pca_t_I)[2]=paste0(vector[i],".pca_t")
        rm(result_pca_t)
      }
      if("pca_wilcoxon" %in% enrichment_method){
        print("pca_wilcoxon will start.")
        result=Get_wilcox_test(geneset=geneset,group=group)
        result_pca_wilcoxon=result[,c("geneset_id","P.Value","logFC")]
        result_pca_wilcoxon$two=TRUE
        result_pca_wilcoxon$invert_A=ifelse(result_pca_wilcoxon$logFC<=0,TRUE,FALSE)
        result_pca_wilcoxon$invert_I=ifelse(result_pca_wilcoxon$logFC>0,TRUE,FALSE)
        result_pca_wilcoxon$P.Value.A=two2one(result_pca_wilcoxon$P.Value, two = as.logical(result_pca_wilcoxon$two), invert = as.logical(result_pca_wilcoxon$invert_A))
        result_pca_wilcoxon$P.Value.I=two2one(result_pca_wilcoxon$P.Value, two = as.logical(result_pca_wilcoxon$two), invert = as.logical(result_pca_wilcoxon$invert_I))
        result_pca_wilcoxon_A=result_pca_wilcoxon[,c("geneset_id","P.Value.A")]
        result_pca_wilcoxon_I=result_pca_wilcoxon[,c("geneset_id","P.Value.I")]
        colnames(result_pca_wilcoxon_A)[2]=paste0(vector[i],".pca_wilcoxon")
        colnames(result_pca_wilcoxon_I)[2]=paste0(vector[i],".pca_wilcoxon")
        rm(result_pca_wilcoxon)
      }
      if("pca_kruskal" %in% enrichment_method){
        print("pca_kruskal will start.")
        result=Get_kruskal_test(geneset=geneset,group=group)
        result_pca_kruskal=result[,c("geneset_id","P.Value","logFC")]
        result_pca_kruskal$two=TRUE
        result_pca_kruskal$invert_A=ifelse(result_pca_kruskal$logFC<=0,TRUE,FALSE)
        result_pca_kruskal$invert_I=ifelse(result_pca_kruskal$logFC>0,TRUE,FALSE)
        result_pca_kruskal$P.Value.A=two2one(result_pca_kruskal$P.Value, two = as.logical(result_pca_kruskal$two), invert = as.logical(result_pca_kruskal$invert_A))
        result_pca_kruskal$P.Value.I=two2one(result_pca_kruskal$P.Value, two = as.logical(result_pca_kruskal$two), invert = as.logical(result_pca_kruskal$invert_I))
        result_pca_kruskal_A=result_pca_kruskal[,c("geneset_id","P.Value.A")]
        result_pca_kruskal_I=result_pca_kruskal[,c("geneset_id","P.Value.I")]
        colnames(result_pca_kruskal_A)[2]=paste0(vector[i],".pca_kruskal")
        colnames(result_pca_kruskal_I)[2]=paste0(vector[i],".pca_kruskal")
        rm(result_pca_kruskal)
      }
      if("pca_permutation" %in% enrichment_method){
        print("pca_permutation will start.")
        result=Get_permutation_test(geneset=geneset,group=group)
        result_pca_permutation=result[,c("geneset_id","P.Value","logFC")]
        result_pca_permutation$two=TRUE
        result_pca_permutation$invert_A=ifelse(result_pca_permutation$logFC<=0,TRUE,FALSE)
        result_pca_permutation$invert_I=ifelse(result_pca_permutation$logFC>0,TRUE,FALSE)
        result_pca_permutation=na.omit(result_pca_permutation)
        result_pca_permutation$P.Value.A=two2one(result_pca_permutation$P.Value, two = as.logical(result_pca_permutation$two), invert = as.logical(result_pca_permutation$invert_A))
        result_pca_permutation$P.Value.I=two2one(result_pca_permutation$P.Value, two = as.logical(result_pca_permutation$two), invert = as.logical(result_pca_permutation$invert_I))
        result_pca_permutation_A=result_pca_permutation[,c("geneset_id","P.Value.A")]
        result_pca_permutation_I=result_pca_permutation[,c("geneset_id","P.Value.I")]
        colnames(result_pca_permutation_A)[2]=paste0(vector[i],".pca_permutation")
        colnames(result_pca_permutation_I)[2]=paste0(vector[i],".pca_permutation")
        rm(result_pca_permutation)
      }
      if("pca_limma" %in% enrichment_method){
        print("pca_limma will start.")
        result=Get_limma(geneset=geneset,group=group)
        result_pca_limma=result[,c("geneset_id","P.Value","logFC")]
        result_pca_limma$two=TRUE
        result_pca_limma$invert_A=ifelse(result_pca_limma$logFC<=0,TRUE,FALSE)
        result_pca_limma$invert_I=ifelse(result_pca_limma$logFC>0,TRUE,FALSE)
        result_pca_limma$P.Value.A=two2one(result_pca_limma$P.Value, two = as.logical(result_pca_limma$two), invert = as.logical(result_pca_limma$invert_A))
        result_pca_limma$P.Value.I=two2one(result_pca_limma$P.Value, two = as.logical(result_pca_limma$two), invert = as.logical(result_pca_limma$invert_I))
        result_pca_limma_A=result_pca_limma[,c("geneset_id","P.Value.A")]
        result_pca_limma_I=result_pca_limma[,c("geneset_id","P.Value.I")]
        colnames(result_pca_limma_A)[2]=paste0(vector[i],".pca_limma")
        colnames(result_pca_limma_I)[2]=paste0(vector[i],".pca_limma")
        rm(result_pca_limma)
      }
      if("pca_anova" %in% enrichment_method){
        print("pca_anova will start.")
        result=Get_anova(geneset=geneset,group=group)
        result_pca_anova=result[,c("geneset_id","P.Value","logFC")]
        result_pca_anova$two=TRUE
        result_pca_anova$invert_A=ifelse(result_pca_anova$logFC<=0,TRUE,FALSE)
        result_pca_anova$invert_I=ifelse(result_pca_anova$logFC>0,TRUE,FALSE)
        result_pca_anova$P.Value.A=two2one(result_pca_anova$P.Value, two = as.logical(result_pca_anova$two), invert = as.logical(result_pca_anova$invert_A))
        result_pca_anova$P.Value.I=two2one(result_pca_anova$P.Value, two = as.logical(result_pca_anova$two), invert = as.logical(result_pca_anova$invert_I))
        result_pca_anova_A=result_pca_anova[,c("geneset_id","P.Value.A")]
        result_pca_anova_I=result_pca_anova[,c("geneset_id","P.Value.I")]
        colnames(result_pca_anova_A)[2]=paste0(vector[i],".pca_anova")
        colnames(result_pca_anova_I)[2]=paste0(vector[i],".pca_anova")
        rm(result_pca_anova)
      }
    } else {
      print("Not using pca.")
    }
    ##########################################aucell
    methods_to_check=c("aucell_t", "aucell_limma", "aucell_anova", "aucell_wilcoxon", "aucell_permutation", "aucell_kruskal")
    if (any(methods_to_check %in% enrichment_method)){
      geneset=AUCell_run(exprMat=as.matrix(exp),geneSets=geneSets)
      geneset=geneset@assays@data@listData[["AUC"]]
      geneset=as.data.frame(geneset)
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("aucell_t" %in% enrichment_method){
        print("aucell_t will start.")
        result=Get_T_test(geneset=geneset,group=group)
        result_aucell_t=result[,c("geneset_id","P.Value","logFC")]
        result_aucell_t$two=TRUE
        result_aucell_t$invert_A=ifelse(result_aucell_t$logFC<=0,TRUE,FALSE)
        result_aucell_t$invert_I=ifelse(result_aucell_t$logFC>0,TRUE,FALSE)
        result_aucell_t$P.Value.A=two2one(result_aucell_t$P.Value, two = as.logical(result_aucell_t$two), invert = as.logical(result_aucell_t$invert_A))
        result_aucell_t$P.Value.I=two2one(result_aucell_t$P.Value, two = as.logical(result_aucell_t$two), invert = as.logical(result_aucell_t$invert_I))
        result_aucell_t_A=result_aucell_t[,c("geneset_id","P.Value.A")]
        result_aucell_t_I=result_aucell_t[,c("geneset_id","P.Value.I")]
        colnames(result_aucell_t_A)[2]=paste0(vector[i],".aucell_t")
        colnames(result_aucell_t_I)[2]=paste0(vector[i],".aucell_t")
        rm(result_aucell_t)
      }
      if("aucell_wilcoxon" %in% enrichment_method){
        print("aucell_wilcoxon will start.")
        result=Get_wilcox_test(geneset=geneset,group=group)
        result_aucell_wilcoxon=result[,c("geneset_id","P.Value","logFC")]
        result_aucell_wilcoxon$two=TRUE
        result_aucell_wilcoxon$invert_A=ifelse(result_aucell_wilcoxon$logFC<=0,TRUE,FALSE)
        result_aucell_wilcoxon$invert_I=ifelse(result_aucell_wilcoxon$logFC>0,TRUE,FALSE)
        result_aucell_wilcoxon$P.Value.A=two2one(result_aucell_wilcoxon$P.Value, two = as.logical(result_aucell_wilcoxon$two), invert = as.logical(result_aucell_wilcoxon$invert_A))
        result_aucell_wilcoxon$P.Value.I=two2one(result_aucell_wilcoxon$P.Value, two = as.logical(result_aucell_wilcoxon$two), invert = as.logical(result_aucell_wilcoxon$invert_I))
        result_aucell_wilcoxon_A=result_aucell_wilcoxon[,c("geneset_id","P.Value.A")]
        result_aucell_wilcoxon_I=result_aucell_wilcoxon[,c("geneset_id","P.Value.I")]
        colnames(result_aucell_wilcoxon_A)[2]=paste0(vector[i],".aucell_wilcoxon")
        colnames(result_aucell_wilcoxon_I)[2]=paste0(vector[i],".aucell_wilcoxon")
        rm(result_aucell_wilcoxon)
      }
      if("aucell_kruskal" %in% enrichment_method){
        print("aucell_kruskal will start.")
        result=Get_kruskal_test(geneset=geneset,group=group)
        result_aucell_kruskal=result[,c("geneset_id","P.Value","logFC")]
        result_aucell_kruskal$two=TRUE
        result_aucell_kruskal$invert_A=ifelse(result_aucell_kruskal$logFC<=0,TRUE,FALSE)
        result_aucell_kruskal$invert_I=ifelse(result_aucell_kruskal$logFC>0,TRUE,FALSE)
        result_aucell_kruskal$P.Value.A=two2one(result_aucell_kruskal$P.Value, two = as.logical(result_aucell_kruskal$two), invert = as.logical(result_aucell_kruskal$invert_A))
        result_aucell_kruskal$P.Value.I=two2one(result_aucell_kruskal$P.Value, two = as.logical(result_aucell_kruskal$two), invert = as.logical(result_aucell_kruskal$invert_I))
        result_aucell_kruskal_A=result_aucell_kruskal[,c("geneset_id","P.Value.A")]
        result_aucell_kruskal_I=result_aucell_kruskal[,c("geneset_id","P.Value.I")]
        colnames(result_aucell_kruskal_A)[2]=paste0(vector[i],".aucell_kruskal")
        colnames(result_aucell_kruskal_I)[2]=paste0(vector[i],".aucell_kruskal")
        rm(result_aucell_kruskal)
      }
      if("aucell_permutation" %in% enrichment_method){
        print("aucell_permutation will start.")
        result=Get_permutation_test(geneset=geneset,group=group)
        result_aucell_permutation=result[,c("geneset_id","P.Value","logFC")]
        result_aucell_permutation$two=TRUE
        result_aucell_permutation$invert_A=ifelse(result_aucell_permutation$logFC<=0,TRUE,FALSE)
        result_aucell_permutation$invert_I=ifelse(result_aucell_permutation$logFC>0,TRUE,FALSE)
        result_aucell_permutation=na.omit(result_aucell_permutation)
        result_aucell_permutation$P.Value.A=two2one(result_aucell_permutation$P.Value, two = as.logical(result_aucell_permutation$two), invert = as.logical(result_aucell_permutation$invert_A))
        result_aucell_permutation$P.Value.I=two2one(result_aucell_permutation$P.Value, two = as.logical(result_aucell_permutation$two), invert = as.logical(result_aucell_permutation$invert_I))
        result_aucell_permutation_A=result_aucell_permutation[,c("geneset_id","P.Value.A")]
        result_aucell_permutation_I=result_aucell_permutation[,c("geneset_id","P.Value.I")]
        colnames(result_aucell_permutation_A)[2]=paste0(vector[i],".aucell_permutation")
        colnames(result_aucell_permutation_I)[2]=paste0(vector[i],".aucell_permutation")
        rm(result_aucell_permutation)
      }
      if("aucell_limma" %in% enrichment_method){
        print("aucell_limma will start.")
        result=Get_limma(geneset=geneset,group=group)
        result_aucell_limma=result[,c("geneset_id","P.Value","logFC")]
        result_aucell_limma$two=TRUE
        result_aucell_limma$invert_A=ifelse(result_aucell_limma$logFC<=0,TRUE,FALSE)
        result_aucell_limma$invert_I=ifelse(result_aucell_limma$logFC>0,TRUE,FALSE)
        result_aucell_limma$P.Value.A=two2one(result_aucell_limma$P.Value, two = as.logical(result_aucell_limma$two), invert = as.logical(result_aucell_limma$invert_A))
        result_aucell_limma$P.Value.I=two2one(result_aucell_limma$P.Value, two = as.logical(result_aucell_limma$two), invert = as.logical(result_aucell_limma$invert_I))
        result_aucell_limma_A=result_aucell_limma[,c("geneset_id","P.Value.A")]
        result_aucell_limma_I=result_aucell_limma[,c("geneset_id","P.Value.I")]
        colnames(result_aucell_limma_A)[2]=paste0(vector[i],".aucell_limma")
        colnames(result_aucell_limma_I)[2]=paste0(vector[i],".aucell_limma")
        rm(result_aucell_limma)
      }
      if("aucell_anova" %in% enrichment_method){
        print("aucell_anova will start.")
        result=Get_anova(geneset=geneset,group=group)
        result_aucell_anova=result[,c("geneset_id","P.Value","logFC")]
        result_aucell_anova$two=TRUE
        result_aucell_anova$invert_A=ifelse(result_aucell_anova$logFC<=0,TRUE,FALSE)
        result_aucell_anova$invert_I=ifelse(result_aucell_anova$logFC>0,TRUE,FALSE)
        result_aucell_anova$P.Value.A=two2one(result_aucell_anova$P.Value, two = as.logical(result_aucell_anova$two), invert = as.logical(result_aucell_anova$invert_A))
        result_aucell_anova$P.Value.I=two2one(result_aucell_anova$P.Value, two = as.logical(result_aucell_anova$two), invert = as.logical(result_aucell_anova$invert_I))
        result_aucell_anova_A=result_aucell_anova[,c("geneset_id","P.Value.A")]
        result_aucell_anova_I=result_aucell_anova[,c("geneset_id","P.Value.I")]
        colnames(result_aucell_anova_A)[2]=paste0(vector[i],".aucell_anova")
        colnames(result_aucell_anova_I)[2]=paste0(vector[i],".aucell_anova")
        rm(result_aucell_anova)
      }
    } else {
      print("Not using aucell.")
    }
    ##########################################ucell
    methods_to_check=c("ucell_t", "ucell_limma", "ucell_anova", "ucell_wilcoxon", "ucell_permutation", "ucell_kruskal")
    if (any(methods_to_check %in% enrichment_method)){
      geneSets_list=lapply(geneSets, function(x){ x@geneIds })
      names(geneSets_list)=names(geneSets)
      geneset=ScoreSignatures_UCell(as.matrix(exp), geneSets_list)
      geneset=as.data.frame(t(geneset))
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("ucell_t" %in% enrichment_method){
        print("ucell_t will start.")
        result=Get_T_test(geneset=geneset,group=group)
        result_ucell_t=result[,c("geneset_id","P.Value","logFC")]
        result_ucell_t$geneset_id=gsub("_UCell$", "", result_ucell_t$geneset_id)
        result_ucell_t$two=TRUE
        result_ucell_t$invert_A=ifelse(result_ucell_t$logFC<=0,TRUE,FALSE)
        result_ucell_t$invert_I=ifelse(result_ucell_t$logFC>0,TRUE,FALSE)
        result_ucell_t$P.Value.A=two2one(result_ucell_t$P.Value, two = as.logical(result_ucell_t$two), invert = as.logical(result_ucell_t$invert_A))
        result_ucell_t$P.Value.I=two2one(result_ucell_t$P.Value, two = as.logical(result_ucell_t$two), invert = as.logical(result_ucell_t$invert_I))
        result_ucell_t_A=result_ucell_t[,c("geneset_id","P.Value.A")]
        result_ucell_t_I=result_ucell_t[,c("geneset_id","P.Value.I")]
        colnames(result_ucell_t_A)[2]=paste0(vector[i],".ucell_t")
        colnames(result_ucell_t_I)[2]=paste0(vector[i],".ucell_t")
        rm(result_ucell_t)
      }
      if("ucell_wilcoxon" %in% enrichment_method){
        print("ucell_wilcoxon will start.")
        result=Get_wilcox_test(geneset=geneset,group=group)
        result_ucell_wilcoxon=result[,c("geneset_id","P.Value","logFC")]
        result_ucell_wilcoxon$geneset_id=gsub("_UCell$", "", result_ucell_wilcoxon$geneset_id)
        result_ucell_wilcoxon$two=TRUE
        result_ucell_wilcoxon$invert_A=ifelse(result_ucell_wilcoxon$logFC<=0,TRUE,FALSE)
        result_ucell_wilcoxon$invert_I=ifelse(result_ucell_wilcoxon$logFC>0,TRUE,FALSE)
        result_ucell_wilcoxon$P.Value.A=two2one(result_ucell_wilcoxon$P.Value, two = as.logical(result_ucell_wilcoxon$two), invert = as.logical(result_ucell_wilcoxon$invert_A))
        result_ucell_wilcoxon$P.Value.I=two2one(result_ucell_wilcoxon$P.Value, two = as.logical(result_ucell_wilcoxon$two), invert = as.logical(result_ucell_wilcoxon$invert_I))
        result_ucell_wilcoxon_A=result_ucell_wilcoxon[,c("geneset_id","P.Value.A")]
        result_ucell_wilcoxon_I=result_ucell_wilcoxon[,c("geneset_id","P.Value.I")]
        colnames(result_ucell_wilcoxon_A)[2]=paste0(vector[i],".ucell_wilcoxon")
        colnames(result_ucell_wilcoxon_I)[2]=paste0(vector[i],".ucell_wilcoxon")
        rm(result_ucell_wilcoxon)
      }
      if("ucell_kruskal" %in% enrichment_method){
        print("ucell_kruskal will start.")
        result=Get_kruskal_test(geneset=geneset,group=group)
        result_ucell_kruskal=result[,c("geneset_id","P.Value","logFC")]
        result_ucell_kruskal$geneset_id=gsub("_UCell$", "", result_ucell_kruskal$geneset_id)
        result_ucell_kruskal$two=TRUE
        result_ucell_kruskal$invert_A=ifelse(result_ucell_kruskal$logFC<=0,TRUE,FALSE)
        result_ucell_kruskal$invert_I=ifelse(result_ucell_kruskal$logFC>0,TRUE,FALSE)
        result_ucell_kruskal$P.Value.A=two2one(result_ucell_kruskal$P.Value, two = as.logical(result_ucell_kruskal$two), invert = as.logical(result_ucell_kruskal$invert_A))
        result_ucell_kruskal$P.Value.I=two2one(result_ucell_kruskal$P.Value, two = as.logical(result_ucell_kruskal$two), invert = as.logical(result_ucell_kruskal$invert_I))
        result_ucell_kruskal_A=result_ucell_kruskal[,c("geneset_id","P.Value.A")]
        result_ucell_kruskal_I=result_ucell_kruskal[,c("geneset_id","P.Value.I")]
        colnames(result_ucell_kruskal_A)[2]=paste0(vector[i],".ucell_kruskal")
        colnames(result_ucell_kruskal_I)[2]=paste0(vector[i],".ucell_kruskal")
        rm(result_ucell_kruskal)
      }
      if("ucell_permutation" %in% enrichment_method){
        print("ucell_permutation will start.")
        result=Get_permutation_test(geneset=geneset,group=group)
        result_ucell_permutation=result[,c("geneset_id","P.Value","logFC")]
        result_ucell_permutation$geneset_id=gsub("_UCell$", "", result_ucell_permutation$geneset_id)
        result_ucell_permutation$two=TRUE
        result_ucell_permutation$invert_A=ifelse(result_ucell_permutation$logFC<=0,TRUE,FALSE)
        result_ucell_permutation$invert_I=ifelse(result_ucell_permutation$logFC>0,TRUE,FALSE)
        result_ucell_permutation=na.omit(result_ucell_permutation)
        result_ucell_permutation$P.Value.A=two2one(result_ucell_permutation$P.Value, two = as.logical(result_ucell_permutation$two), invert = as.logical(result_ucell_permutation$invert_A))
        result_ucell_permutation$P.Value.I=two2one(result_ucell_permutation$P.Value, two = as.logical(result_ucell_permutation$two), invert = as.logical(result_ucell_permutation$invert_I))
        result_ucell_permutation_A=result_ucell_permutation[,c("geneset_id","P.Value.A")]
        result_ucell_permutation_I=result_ucell_permutation[,c("geneset_id","P.Value.I")]
        colnames(result_ucell_permutation_A)[2]=paste0(vector[i],".ucell_permutation")
        colnames(result_ucell_permutation_I)[2]=paste0(vector[i],".ucell_permutation")
        rm(result_ucell_permutation)
      }
      if("ucell_limma" %in% enrichment_method){
        print("ucell_limma will start.")
        result=Get_limma(geneset=geneset,group=group)
        result_ucell_limma=result[,c("geneset_id","P.Value","logFC")]
        result_ucell_limma$geneset_id=gsub("_UCell$", "", result_ucell_limma$geneset_id)
        result_ucell_limma$two=TRUE
        result_ucell_limma$invert_A=ifelse(result_ucell_limma$logFC<=0,TRUE,FALSE)
        result_ucell_limma$invert_I=ifelse(result_ucell_limma$logFC>0,TRUE,FALSE)
        result_ucell_limma$P.Value.A=two2one(result_ucell_limma$P.Value, two = as.logical(result_ucell_limma$two), invert = as.logical(result_ucell_limma$invert_A))
        result_ucell_limma$P.Value.I=two2one(result_ucell_limma$P.Value, two = as.logical(result_ucell_limma$two), invert = as.logical(result_ucell_limma$invert_I))
        result_ucell_limma_A=result_ucell_limma[,c("geneset_id","P.Value.A")]
        result_ucell_limma_I=result_ucell_limma[,c("geneset_id","P.Value.I")]
        colnames(result_ucell_limma_A)[2]=paste0(vector[i],".ucell_limma")
        colnames(result_ucell_limma_I)[2]=paste0(vector[i],".ucell_limma")
        rm(result_ucell_limma)
      }
      if("ucell_anova" %in% enrichment_method){
        print("ucell_anova will start.")
        result=Get_anova(geneset=geneset,group=group)
        result_ucell_anova=result[,c("geneset_id","P.Value","logFC")]
        result_ucell_anova$geneset_id=gsub("_UCell$", "", result_ucell_anova$geneset_id)
        result_ucell_anova$two=TRUE
        result_ucell_anova$invert_A=ifelse(result_ucell_anova$logFC<=0,TRUE,FALSE)
        result_ucell_anova$invert_I=ifelse(result_ucell_anova$logFC>0,TRUE,FALSE)
        result_ucell_anova$P.Value.A=two2one(result_ucell_anova$P.Value, two = as.logical(result_ucell_anova$two), invert = as.logical(result_ucell_anova$invert_A))
        result_ucell_anova$P.Value.I=two2one(result_ucell_anova$P.Value, two = as.logical(result_ucell_anova$two), invert = as.logical(result_ucell_anova$invert_I))
        result_ucell_anova_A=result_ucell_anova[,c("geneset_id","P.Value.A")]
        result_ucell_anova_I=result_ucell_anova[,c("geneset_id","P.Value.I")]
        colnames(result_ucell_anova_A)[2]=paste0(vector[i],".ucell_anova")
        colnames(result_ucell_anova_I)[2]=paste0(vector[i],".ucell_anova")
        rm(result_ucell_anova)
      }
    } else {
      print("Not using ucell.")
    }
    ##########################################singscore
    methods_to_check=c("singscore_t", "singscore_limma", "singscore_anova", "singscore_wilcoxon", "singscore_permutation", "singscore_kruskal")
    if (any(methods_to_check %in% enrichment_method)){
      geneSets_list=lapply(geneSets, function(x){ x@geneIds })
      names(geneSets_list)=names(geneSets)
      eranks = rankGenes(as.matrix(exp))
      geneset=data.frame()
      for(eps in 1:length(geneSets_list)){
        epi_sig=GeneSet(geneSets_list[[eps]], setName = names(geneSets_list)[eps], geneIdType = SymbolIdentifier())
        epi_score=simpleScore(eranks, epi_sig)
        epi_score=as.data.frame(t(epi_score[,-2,drop=F]))
        rownames(epi_score)=names(geneSets_list)[eps]
        geneset=rbind(geneset,epi_score)
      }
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("singscore_t" %in% enrichment_method){
        print("singscore_t will start.")
        result=Get_T_test(geneset=geneset,group=group)
        result_singscore_t=result[,c("geneset_id","P.Value","logFC")]
        result_singscore_t$two=TRUE
        result_singscore_t$invert_A=ifelse(result_singscore_t$logFC<=0,TRUE,FALSE)
        result_singscore_t$invert_I=ifelse(result_singscore_t$logFC>0,TRUE,FALSE)
        result_singscore_t$P.Value.A=two2one(result_singscore_t$P.Value, two = as.logical(result_singscore_t$two), invert = as.logical(result_singscore_t$invert_A))
        result_singscore_t$P.Value.I=two2one(result_singscore_t$P.Value, two = as.logical(result_singscore_t$two), invert = as.logical(result_singscore_t$invert_I))
        result_singscore_t_A=result_singscore_t[,c("geneset_id","P.Value.A")]
        result_singscore_t_I=result_singscore_t[,c("geneset_id","P.Value.I")]
        colnames(result_singscore_t_A)[2]=paste0(vector[i],".singscore_t")
        colnames(result_singscore_t_I)[2]=paste0(vector[i],".singscore_t")
        rm(result_singscore_t)
      }
      if("singscore_wilcoxon" %in% enrichment_method){
        print("singscore_wilcoxon will start.")
        result=Get_wilcox_test(geneset=geneset,group=group)
        result_singscore_wilcoxon=result[,c("geneset_id","P.Value","logFC")]
        result_singscore_wilcoxon$two=TRUE
        result_singscore_wilcoxon$invert_A=ifelse(result_singscore_wilcoxon$logFC<=0,TRUE,FALSE)
        result_singscore_wilcoxon$invert_I=ifelse(result_singscore_wilcoxon$logFC>0,TRUE,FALSE)
        result_singscore_wilcoxon$P.Value.A=two2one(result_singscore_wilcoxon$P.Value, two = as.logical(result_singscore_wilcoxon$two), invert = as.logical(result_singscore_wilcoxon$invert_A))
        result_singscore_wilcoxon$P.Value.I=two2one(result_singscore_wilcoxon$P.Value, two = as.logical(result_singscore_wilcoxon$two), invert = as.logical(result_singscore_wilcoxon$invert_I))
        result_singscore_wilcoxon_A=result_singscore_wilcoxon[,c("geneset_id","P.Value.A")]
        result_singscore_wilcoxon_I=result_singscore_wilcoxon[,c("geneset_id","P.Value.I")]
        colnames(result_singscore_wilcoxon_A)[2]=paste0(vector[i],".singscore_wilcoxon")
        colnames(result_singscore_wilcoxon_I)[2]=paste0(vector[i],".singscore_wilcoxon")
        rm(result_singscore_wilcoxon)
      }
      if("singscore_kruskal" %in% enrichment_method){
        print("singscore_kruskal will start.")
        result=Get_kruskal_test(geneset=geneset,group=group)
        result_singscore_kruskal=result[,c("geneset_id","P.Value","logFC")]
        result_singscore_kruskal$two=TRUE
        result_singscore_kruskal$invert_A=ifelse(result_singscore_kruskal$logFC<=0,TRUE,FALSE)
        result_singscore_kruskal$invert_I=ifelse(result_singscore_kruskal$logFC>0,TRUE,FALSE)
        result_singscore_kruskal$P.Value.A=two2one(result_singscore_kruskal$P.Value, two = as.logical(result_singscore_kruskal$two), invert = as.logical(result_singscore_kruskal$invert_A))
        result_singscore_kruskal$P.Value.I=two2one(result_singscore_kruskal$P.Value, two = as.logical(result_singscore_kruskal$two), invert = as.logical(result_singscore_kruskal$invert_I))
        result_singscore_kruskal_A=result_singscore_kruskal[,c("geneset_id","P.Value.A")]
        result_singscore_kruskal_I=result_singscore_kruskal[,c("geneset_id","P.Value.I")]
        colnames(result_singscore_kruskal_A)[2]=paste0(vector[i],".singscore_kruskal")
        colnames(result_singscore_kruskal_I)[2]=paste0(vector[i],".singscore_kruskal")
        rm(result_singscore_kruskal)
      }
      if("singscore_permutation" %in% enrichment_method){
        print("singscore_permutation will start.")
        result=Get_permutation_test(geneset=geneset,group=group)
        result_singscore_permutation=result[,c("geneset_id","P.Value","logFC")]
        result_singscore_permutation$two=TRUE
        result_singscore_permutation$invert_A=ifelse(result_singscore_permutation$logFC<=0,TRUE,FALSE)
        result_singscore_permutation$invert_I=ifelse(result_singscore_permutation$logFC>0,TRUE,FALSE)
        result_singscore_permutation=na.omit(result_singscore_permutation)
        result_singscore_permutation$P.Value.A=two2one(result_singscore_permutation$P.Value, two = as.logical(result_singscore_permutation$two), invert = as.logical(result_singscore_permutation$invert_A))
        result_singscore_permutation$P.Value.I=two2one(result_singscore_permutation$P.Value, two = as.logical(result_singscore_permutation$two), invert = as.logical(result_singscore_permutation$invert_I))
        result_singscore_permutation_A=result_singscore_permutation[,c("geneset_id","P.Value.A")]
        result_singscore_permutation_I=result_singscore_permutation[,c("geneset_id","P.Value.I")]
        colnames(result_singscore_permutation_A)[2]=paste0(vector[i],".singscore_permutation")
        colnames(result_singscore_permutation_I)[2]=paste0(vector[i],".singscore_permutation")
        rm(result_singscore_permutation)
      }
      if("singscore_limma" %in% enrichment_method){
        print("singscore_limma will start.")
        result=Get_limma(geneset=geneset,group=group)
        result_singscore_limma=result[,c("geneset_id","P.Value","logFC")]
        result_singscore_limma$two=TRUE
        result_singscore_limma$invert_A=ifelse(result_singscore_limma$logFC<=0,TRUE,FALSE)
        result_singscore_limma$invert_I=ifelse(result_singscore_limma$logFC>0,TRUE,FALSE)
        result_singscore_limma$P.Value.A=two2one(result_singscore_limma$P.Value, two = as.logical(result_singscore_limma$two), invert = as.logical(result_singscore_limma$invert_A))
        result_singscore_limma$P.Value.I=two2one(result_singscore_limma$P.Value, two = as.logical(result_singscore_limma$two), invert = as.logical(result_singscore_limma$invert_I))
        result_singscore_limma_A=result_singscore_limma[,c("geneset_id","P.Value.A")]
        result_singscore_limma_I=result_singscore_limma[,c("geneset_id","P.Value.I")]
        colnames(result_singscore_limma_A)[2]=paste0(vector[i],".singscore_limma")
        colnames(result_singscore_limma_I)[2]=paste0(vector[i],".singscore_limma")
        rm(result_singscore_limma)
      }
      if("singscore_anova" %in% enrichment_method){
        print("singscore_anova will start.")
        result=Get_anova(geneset=geneset,group=group)
        result_singscore_anova=result[,c("geneset_id","P.Value","logFC")]
        result_singscore_anova$two=TRUE
        result_singscore_anova$invert_A=ifelse(result_singscore_anova$logFC<=0,TRUE,FALSE)
        result_singscore_anova$invert_I=ifelse(result_singscore_anova$logFC>0,TRUE,FALSE)
        result_singscore_anova$P.Value.A=two2one(result_singscore_anova$P.Value, two = as.logical(result_singscore_anova$two), invert = as.logical(result_singscore_anova$invert_A))
        result_singscore_anova$P.Value.I=two2one(result_singscore_anova$P.Value, two = as.logical(result_singscore_anova$two), invert = as.logical(result_singscore_anova$invert_I))
        result_singscore_anova_A=result_singscore_anova[,c("geneset_id","P.Value.A")]
        result_singscore_anova_I=result_singscore_anova[,c("geneset_id","P.Value.I")]
        colnames(result_singscore_anova_A)[2]=paste0(vector[i],".singscore_anova")
        colnames(result_singscore_anova_I)[2]=paste0(vector[i],".singscore_anova")
        rm(result_singscore_anova)
      }
    } else {
      print("Not using singscore.")
    }
    ##########################################median
    methods_to_check=c("median_t", "median_limma", "median_anova", "median_wilcoxon", "median_permutation", "median_kruskal")
    if (any(methods_to_check %in% enrichment_method)){
      geneset=Median_evaluation(exp=exp,geneSets=geneSets)
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("median_t" %in% enrichment_method){
        print("median_t will start.")
        result=Get_T_test(geneset=geneset,group=group)
        result_median_t=result[,c("geneset_id","P.Value","logFC")]
        result_median_t$two=TRUE
        result_median_t$invert_A=ifelse(result_median_t$logFC<=0,TRUE,FALSE)
        result_median_t$invert_I=ifelse(result_median_t$logFC>0,TRUE,FALSE)
        result_median_t$P.Value.A=two2one(result_median_t$P.Value, two = as.logical(result_median_t$two), invert = as.logical(result_median_t$invert_A))
        result_median_t$P.Value.I=two2one(result_median_t$P.Value, two = as.logical(result_median_t$two), invert = as.logical(result_median_t$invert_I))
        result_median_t_A=result_median_t[,c("geneset_id","P.Value.A")]
        result_median_t_I=result_median_t[,c("geneset_id","P.Value.I")]
        colnames(result_median_t_A)[2]=paste0(vector[i],".median_t")
        colnames(result_median_t_I)[2]=paste0(vector[i],".median_t")
        rm(result_median_t)
      }
      if("median_wilcoxon" %in% enrichment_method){
        print("median_wilcoxon will start.")
        result=Get_wilcox_test(geneset=geneset,group=group)
        result_median_wilcoxon=result[,c("geneset_id","P.Value","logFC")]
        result_median_wilcoxon$two=TRUE
        result_median_wilcoxon$invert_A=ifelse(result_median_wilcoxon$logFC<=0,TRUE,FALSE)
        result_median_wilcoxon$invert_I=ifelse(result_median_wilcoxon$logFC>0,TRUE,FALSE)
        result_median_wilcoxon$P.Value.A=two2one(result_median_wilcoxon$P.Value, two = as.logical(result_median_wilcoxon$two), invert = as.logical(result_median_wilcoxon$invert_A))
        result_median_wilcoxon$P.Value.I=two2one(result_median_wilcoxon$P.Value, two = as.logical(result_median_wilcoxon$two), invert = as.logical(result_median_wilcoxon$invert_I))
        result_median_wilcoxon_A=result_median_wilcoxon[,c("geneset_id","P.Value.A")]
        result_median_wilcoxon_I=result_median_wilcoxon[,c("geneset_id","P.Value.I")]
        colnames(result_median_wilcoxon_A)[2]=paste0(vector[i],".median_wilcoxon")
        colnames(result_median_wilcoxon_I)[2]=paste0(vector[i],".median_wilcoxon")
        rm(result_median_wilcoxon)
      }
      if("median_kruskal" %in% enrichment_method){
        print("median_kruskal will start.")
        result=Get_kruskal_test(geneset=geneset,group=group)
        result_median_kruskal=result[,c("geneset_id","P.Value","logFC")]
        result_median_kruskal$two=TRUE
        result_median_kruskal$invert_A=ifelse(result_median_kruskal$logFC<=0,TRUE,FALSE)
        result_median_kruskal$invert_I=ifelse(result_median_kruskal$logFC>0,TRUE,FALSE)
        result_median_kruskal$P.Value.A=two2one(result_median_kruskal$P.Value, two = as.logical(result_median_kruskal$two), invert = as.logical(result_median_kruskal$invert_A))
        result_median_kruskal$P.Value.I=two2one(result_median_kruskal$P.Value, two = as.logical(result_median_kruskal$two), invert = as.logical(result_median_kruskal$invert_I))
        result_median_kruskal_A=result_median_kruskal[,c("geneset_id","P.Value.A")]
        result_median_kruskal_I=result_median_kruskal[,c("geneset_id","P.Value.I")]
        colnames(result_median_kruskal_A)[2]=paste0(vector[i],".median_kruskal")
        colnames(result_median_kruskal_I)[2]=paste0(vector[i],".median_kruskal")
        rm(result_median_kruskal)
      }
      if("median_permutation" %in% enrichment_method){
        print("median_permutation will start.")
        result=Get_permutation_test(geneset=geneset,group=group)
        result_median_permutation=result[,c("geneset_id","P.Value","logFC")]
        result_median_permutation$two=TRUE
        result_median_permutation$invert_A=ifelse(result_median_permutation$logFC<=0,TRUE,FALSE)
        result_median_permutation$invert_I=ifelse(result_median_permutation$logFC>0,TRUE,FALSE)
        result_median_permutation=na.omit(result_median_permutation)
        result_median_permutation$P.Value.A=two2one(result_median_permutation$P.Value, two = as.logical(result_median_permutation$two), invert = as.logical(result_median_permutation$invert_A))
        result_median_permutation$P.Value.I=two2one(result_median_permutation$P.Value, two = as.logical(result_median_permutation$two), invert = as.logical(result_median_permutation$invert_I))
        result_median_permutation_A=result_median_permutation[,c("geneset_id","P.Value.A")]
        result_median_permutation_I=result_median_permutation[,c("geneset_id","P.Value.I")]
        colnames(result_median_permutation_A)[2]=paste0(vector[i],".median_permutation")
        colnames(result_median_permutation_I)[2]=paste0(vector[i],".median_permutation")
        rm(result_median_permutation)
      }
      if("median_limma" %in% enrichment_method){
        print("median_limma will start.")
        result=Get_limma(geneset=geneset,group=group)
        result_median_limma=result[,c("geneset_id","P.Value","logFC")]
        result_median_limma$two=TRUE
        result_median_limma$invert_A=ifelse(result_median_limma$logFC<=0,TRUE,FALSE)
        result_median_limma$invert_I=ifelse(result_median_limma$logFC>0,TRUE,FALSE)
        result_median_limma$P.Value.A=two2one(result_median_limma$P.Value, two = as.logical(result_median_limma$two), invert = as.logical(result_median_limma$invert_A))
        result_median_limma$P.Value.I=two2one(result_median_limma$P.Value, two = as.logical(result_median_limma$two), invert = as.logical(result_median_limma$invert_I))
        result_median_limma_A=result_median_limma[,c("geneset_id","P.Value.A")]
        result_median_limma_I=result_median_limma[,c("geneset_id","P.Value.I")]
        colnames(result_median_limma_A)[2]=paste0(vector[i],".median_limma")
        colnames(result_median_limma_I)[2]=paste0(vector[i],".median_limma")
        rm(result_median_limma)
      }
      if("median_anova" %in% enrichment_method){
        print("median_anova will start.")
        result=Get_anova(geneset=geneset,group=group)
        result_median_anova=result[,c("geneset_id","P.Value","logFC")]
        result_median_anova$two=TRUE
        result_median_anova$invert_A=ifelse(result_median_anova$logFC<=0,TRUE,FALSE)
        result_median_anova$invert_I=ifelse(result_median_anova$logFC>0,TRUE,FALSE)
        result_median_anova$P.Value.A=two2one(result_median_anova$P.Value, two = as.logical(result_median_anova$two), invert = as.logical(result_median_anova$invert_A))
        result_median_anova$P.Value.I=two2one(result_median_anova$P.Value, two = as.logical(result_median_anova$two), invert = as.logical(result_median_anova$invert_I))
        result_median_anova_A=result_median_anova[,c("geneset_id","P.Value.A")]
        result_median_anova_I=result_median_anova[,c("geneset_id","P.Value.I")]
        colnames(result_median_anova_A)[2]=paste0(vector[i],".median_anova")
        colnames(result_median_anova_I)[2]=paste0(vector[i],".median_anova")
        rm(result_median_anova)
      }
    } else {
      print("Not using median.")
    }
    ##########################################fgsea/ora
    methods_to_check=c("t_fgsea", "limma_fgsea", "anova_fgsea", "wilcoxon_fgsea", "permutation_fgsea", "kruskal_fgsea",
                       "t_ora","limma_ora","anova_ora","wilcoxon_ora","permutation_ora","kruskal_ora")
    if (any(methods_to_check %in% enrichment_method)){
      if (any(c("t_fgsea","t_ora") %in% enrichment_method)){
        all_diff=Get_T_test(geneset=exp,group=group)
        if("t_fgsea"%in% enrichment_method){
          print("t_fgsea will start.")
          alldiff=all_diff[,c("geneset_id","P.Value","logFC")]
          alldiff$id=-log10(alldiff$P.Value)*alldiff$logFC
          alldiff=alldiff[order(alldiff$id,decreasing = T),]
          id=alldiff$id
          names(id)=alldiff$geneset_id
          canonical_pathways=do.call(rbind, lapply(geneSets, function(gs) {  
            data.frame(term = gs@setName, gene = gs@geneIds, stringsAsFactors = FALSE)  
          }))
          canonical_pathways=canonical_pathways[canonical_pathways$gene != "", ]
          canonical_pathways$term=as.character(canonical_pathways$term)
          canonical_pathways=canonical_pathways %>% split(x = .$gene, f = .$term)
          fgseaRes=fgsea(pathways=canonical_pathways,
                         stats=id[is.finite(id)],
                         minSize=min.sz,
                         maxSize=max.sz,
                         nperm=1000)
          result=data.frame(fgseaRes)
          result_t_fgsea=result[,c("pathway","pval","ES")]
          colnames(result_t_fgsea)=c("geneset_id","P.Value","ES")
          result_t_fgsea$geneset_id=gsub("\\.", "'", result_t_fgsea$geneset_id)
          result_t_fgsea=na.omit(result_t_fgsea)
          result_t_fgsea$two=TRUE
          result_t_fgsea$invert_A=ifelse(result_t_fgsea$ES<0,TRUE,FALSE)
          result_t_fgsea$invert_I=ifelse(result_t_fgsea$ES>0,TRUE,FALSE)
          result_t_fgsea$P.Value.A=two2one(result_t_fgsea$P.Value, two = as.logical(result_t_fgsea$two), invert = as.logical(result_t_fgsea$invert_A))
          result_t_fgsea$P.Value.I=two2one(result_t_fgsea$P.Value, two = as.logical(result_t_fgsea$two), invert = as.logical(result_t_fgsea$invert_I))
          result_t_fgsea_A=result_t_fgsea[,c("geneset_id","P.Value.A")]
          result_t_fgsea_I=result_t_fgsea[,c("geneset_id","P.Value.I")]
          colnames(result_t_fgsea_A)[2]=paste0(vector[i],".t_fgsea")
          colnames(result_t_fgsea_I)[2]=paste0(vector[i],".t_fgsea")
          rm(result_t_fgsea)
        }
        if("t_ora"%in% enrichment_method){
          print("t_ora will start.")
          result_UP=all_diff[all_diff$logFC>=0,]
          result_DN=all_diff[all_diff$logFC<0,]
          result_UP=result_UP[order(result_UP$P.Value, -result_UP$logFC), ] 
          result_DN=result_DN[order(result_DN$P.Value, result_DN$logFC), ]
          upregulated_gene_sets=head(result_UP,min(1000, nrow(result_UP)))
          downregulated_gene_sets=head(result_DN,min(1000, nrow(result_DN)))
          gene_H=upregulated_gene_sets$geneset_id
          gene_L=downregulated_gene_sets$geneset_id
          canonical_pathways=do.call(rbind, lapply(geneSets, function(gs) {  
            data.frame(term = gs@setName, gene = gs@geneIds, stringsAsFactors = FALSE)  
          }))
          canonical_pathways=canonical_pathways[canonical_pathways$gene != "", ]
          enrich_H=enricher(gene_H, TERM2GENE = canonical_pathways,minGSSize=min.sz,maxGSSize=max.sz,
                            pvalueCutoff = 1,qvalueCutoff = 1)
          if (is.null(enrich_H)) {
            result_t_ora_A=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_t_ora_A=enrich_H@result
            result_t_ora_A=result_t_ora_A[,c("Description","pvalue")]
          }
          colnames(result_t_ora_A)=c("geneset_id","P.Value_H")
          result_t_ora_A$geneset_id=gsub("\\.", "'", result_t_ora_A$geneset_id)
          enrich_L=enricher(gene_L, TERM2GENE = canonical_pathways,minGSSize=min.sz,maxGSSize=max.sz,
                            pvalueCutoff = 1,qvalueCutoff = 1)
          if (is.null(enrich_L)) {
            result_t_ora_I=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_t_ora_I=enrich_L@result
            result_t_ora_I=result_t_ora_I[,c("Description","pvalue")]
          }
          colnames(result_t_ora_I)=c("geneset_id","P.Value_L")
          result_t_ora_I$geneset_id=gsub("\\.", "'", result_t_ora_I$geneset_id)
          colnames(result_t_ora_A)[2]=paste0(vector[i],".t_ora")
          colnames(result_t_ora_I)[2]=paste0(vector[i],".t_ora")
        }
      }
      if (any(c("limma_fgsea","limma_ora") %in% enrichment_method)){
        all_diff=Get_limma(geneset=exp,group=group)
        if("limma_fgsea" %in% enrichment_method){
          print("limma_fgsea will start.")
          alldiff=all_diff[,c("geneset_id","P.Value","logFC")]
          alldiff$id=-log10(alldiff$P.Value)*alldiff$logFC
          alldiff=alldiff[order(alldiff$id,decreasing = T),]
          id=alldiff$id
          names(id)=alldiff$geneset_id
          canonical_pathways=do.call(rbind, lapply(geneSets, function(gs) {  
            data.frame(term = gs@setName, gene = gs@geneIds, stringsAsFactors = FALSE)  
          }))
          canonical_pathways=canonical_pathways[canonical_pathways$gene != "", ]
          canonical_pathways$term=as.character(canonical_pathways$term)
          canonical_pathways=canonical_pathways %>% split(x = .$gene, f = .$term)
          fgseaRes=fgsea(pathways=canonical_pathways,
                         stats=id[is.finite(id)],
                         minSize=min.sz,
                         maxSize=max.sz,
                         nperm=1000)
          result=data.frame(fgseaRes)
          result_limma_fgsea=result[,c("pathway","pval","ES")]
          colnames(result_limma_fgsea)=c("geneset_id","P.Value","ES")
          result_limma_fgsea$geneset_id=gsub("\\.", "'", result_limma_fgsea$geneset_id)
          result_limma_fgsea=na.omit(result_limma_fgsea)
          result_limma_fgsea$two=TRUE
          result_limma_fgsea$invert_A=ifelse(result_limma_fgsea$ES<0,TRUE,FALSE)
          result_limma_fgsea$invert_I=ifelse(result_limma_fgsea$ES>0,TRUE,FALSE)
          result_limma_fgsea$P.Value.A=two2one(result_limma_fgsea$P.Value, two = as.logical(result_limma_fgsea$two), invert = as.logical(result_limma_fgsea$invert_A))
          result_limma_fgsea$P.Value.I=two2one(result_limma_fgsea$P.Value, two = as.logical(result_limma_fgsea$two), invert = as.logical(result_limma_fgsea$invert_I))
          result_limma_fgsea_A=result_limma_fgsea[,c("geneset_id","P.Value.A")]
          result_limma_fgsea_I=result_limma_fgsea[,c("geneset_id","P.Value.I")]
          colnames(result_limma_fgsea_A)[2]=paste0(vector[i],".limma_fgsea")
          colnames(result_limma_fgsea_I)[2]=paste0(vector[i],".limma_fgsea")
          rm(result_limma_fgsea)
        }
        if("limma_ora" %in% enrichment_method){
          print("limma_ora will start.")
          result_UP=all_diff[all_diff$logFC>=0,]
          result_DN=all_diff[all_diff$logFC<0,]
          result_UP=result_UP[order(result_UP$P.Value, -result_UP$logFC), ] 
          result_DN=result_DN[order(result_DN$P.Value, result_DN$logFC), ]
          upregulated_gene_sets=head(result_UP,min(1000, nrow(result_UP)))
          downregulated_gene_sets=head(result_DN,min(1000, nrow(result_DN)))
          gene_H=rownames(upregulated_gene_sets)
          gene_L=rownames(downregulated_gene_sets)
          canonical_pathways=do.call(rbind, lapply(geneSets, function(gs) {  
            data.frame(term = gs@setName, gene = gs@geneIds, stringsAsFactors = FALSE)  
          }))
          canonical_pathways=canonical_pathways[canonical_pathways$gene != "", ]
          enrich_H=enricher(gene_H, TERM2GENE=canonical_pathways,minGSSize = min.sz,maxGSSize=max.sz,
                            pvalueCutoff = 1,qvalueCutoff = 1)
          if (is.null(enrich_H)) {
            result_limma_ora_A=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_limma_ora_A=enrich_H@result
            result_limma_ora_A=result_limma_ora_A[,c("Description","pvalue")]
          }
          colnames(result_limma_ora_A)=c("geneset_id","P.Value_H")
          result_limma_ora_A$geneset_id=gsub("\\.", "'", result_limma_ora_A$geneset_id)
          enrich_L=enricher(gene_L, TERM2GENE=canonical_pathways,minGSSize = min.sz,maxGSSize=max.sz,
                            pvalueCutoff=1,qvalueCutoff=1)
          if (is.null(enrich_L)) {
            result_limma_ora_I=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_limma_ora_I=enrich_L@result
            result_limma_ora_I=result_limma_ora_I[,c("Description","pvalue")]
          }
          colnames(result_limma_ora_I)=c("geneset_id","P.Value_L")
          result_limma_ora_I$geneset_id=gsub("\\.", "'", result_limma_ora_I$geneset_id)
          colnames(result_limma_ora_A)[2]=paste0(vector[i],".limma_ora")
          colnames(result_limma_ora_I)[2]=paste0(vector[i],".limma_ora")
        }
      }
      if (any(c("anova_fgsea","anova_ora") %in% enrichment_method)){
        all_diff=Get_anova(geneset=exp,group=group)
        if("anova_fgsea" %in% enrichment_method){
          print("anova_fgsea will start.")
          alldiff=all_diff[,c("geneset_id","P.Value","logFC")]
          alldiff$id=-log10(alldiff$P.Value)*alldiff$logFC
          alldiff=alldiff[order(alldiff$id,decreasing = T),]
          id=alldiff$id
          names(id)=alldiff$geneset_id
          canonical_pathways=do.call(rbind, lapply(geneSets, function(gs) {  
            data.frame(term = gs@setName, gene = gs@geneIds, stringsAsFactors = FALSE)  
          }))
          canonical_pathways=canonical_pathways[canonical_pathways$gene != "", ]
          canonical_pathways$term=as.character(canonical_pathways$term)
          canonical_pathways=canonical_pathways %>% split(x = .$gene, f = .$term)
          fgseaRes=fgsea(pathways=canonical_pathways,
                         stats=id[is.finite(id)],
                         minSize=min.sz,
                         maxSize=max.sz,
                         nperm=1000)
          result=data.frame(fgseaRes)
          result_anova_fgsea=result[,c("pathway","pval","ES")]
          colnames(result_anova_fgsea)=c("geneset_id","P.Value","ES")
          result_anova_fgsea$geneset_id=gsub("\\.", "'", result_anova_fgsea$geneset_id)
          result_anova_fgsea=na.omit(result_anova_fgsea)
          result_anova_fgsea$two=TRUE
          result_anova_fgsea$invert_A=ifelse(result_anova_fgsea$ES<0,TRUE,FALSE)
          result_anova_fgsea$invert_I=ifelse(result_anova_fgsea$ES>0,TRUE,FALSE)
          result_anova_fgsea$P.Value.A=two2one(result_anova_fgsea$P.Value, two = as.logical(result_anova_fgsea$two), invert = as.logical(result_anova_fgsea$invert_A))
          result_anova_fgsea$P.Value.I=two2one(result_anova_fgsea$P.Value, two = as.logical(result_anova_fgsea$two), invert = as.logical(result_anova_fgsea$invert_I))
          result_anova_fgsea_A=result_anova_fgsea[,c("geneset_id","P.Value.A")]
          result_anova_fgsea_I=result_anova_fgsea[,c("geneset_id","P.Value.I")]
          colnames(result_anova_fgsea_A)[2]=paste0(vector[i],".anova_fgsea")
          colnames(result_anova_fgsea_I)[2]=paste0(vector[i],".anova_fgsea")
          rm(result_anova_fgsea)
        }
        if("anova_ora" %in% enrichment_method){
          print("anova_ora will start.")
          result_UP=all_diff[all_diff$logFC>=0,]
          result_DN=all_diff[all_diff$logFC<0,]
          result_UP=result_UP[order(result_UP$P.Value, -result_UP$logFC), ] 
          result_DN=result_DN[order(result_DN$P.Value, result_DN$logFC), ]
          upregulated_gene_sets=head(result_UP,min(1000, nrow(result_UP)))
          downregulated_gene_sets=head(result_DN,min(1000, nrow(result_DN)))
          gene_H=upregulated_gene_sets$geneset_id
          gene_L=downregulated_gene_sets$geneset_id
          canonical_pathways=do.call(rbind, lapply(geneSets, function(gs) {  
            data.frame(term = gs@setName, gene = gs@geneIds, stringsAsFactors = FALSE)  
          }))
          canonical_pathways=canonical_pathways[canonical_pathways$gene != "", ]
          enrich_H=enricher(gene_H, TERM2GENE = canonical_pathways,minGSSize=min.sz,maxGSSize=max.sz,
                            pvalueCutoff = 1,qvalueCutoff = 1)
          if (is.null(enrich_H)) {
            result_anova_ora_A=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_anova_ora_A=enrich_H@result
            result_anova_ora_A=result_anova_ora_A[,c("Description","pvalue")]
          }
          colnames(result_anova_ora_A)=c("geneset_id","P.Value_H")
          result_anova_ora_A$geneset_id=gsub("\\.", "'", result_anova_ora_A$geneset_id)
          enrich_L=enricher(gene_L, TERM2GENE = canonical_pathways,minGSSize=min.sz,maxGSSize=max.sz,
                            pvalueCutoff = 1,qvalueCutoff = 1)
          if (is.null(enrich_L)) {
            result_anova_ora_I=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_anova_ora_I=enrich_L@result
            result_anova_ora_I=result_anova_ora_I[,c("Description","pvalue")]
          }
          colnames(result_anova_ora_I)=c("geneset_id","P.Value_L")
          result_anova_ora_I$geneset_id=gsub("\\.", "'", result_anova_ora_I$geneset_id)
          colnames(result_anova_ora_A)[2]=paste0(vector[i],".anova_ora")
          colnames(result_anova_ora_I)[2]=paste0(vector[i],".anova_ora")
        }
      }
      if (any(c("wilcoxon_fgsea","wilcoxon_ora") %in% enrichment_method)){
        all_diff=Get_wilcox_test(geneset=exp,group=group)
        if("wilcoxon_fgsea"%in% enrichment_method){
          print("wilcoxon_fgsea will start.")
          alldiff=all_diff[,c("geneset_id","P.Value","logFC")]
          alldiff$id=-log10(alldiff$P.Value)*alldiff$logFC
          alldiff=alldiff[order(alldiff$id,decreasing = T),]
          id=alldiff$id
          names(id)=alldiff$geneset_id
          canonical_pathways=do.call(rbind, lapply(geneSets, function(gs) {  
            data.frame(term = gs@setName, gene = gs@geneIds, stringsAsFactors = FALSE)  
          }))
          canonical_pathways=canonical_pathways[canonical_pathways$gene != "", ]
          canonical_pathways$term=as.character(canonical_pathways$term)
          canonical_pathways=canonical_pathways %>% split(x = .$gene, f = .$term)
          fgseaRes=fgsea(pathways=canonical_pathways,
                         stats=id[is.finite(id)],
                         minSize=min.sz,
                         maxSize=max.sz,
                         nperm=1000)
          result=data.frame(fgseaRes)
          result_wilcoxon_fgsea=result[,c("pathway","pval","ES")]
          colnames(result_wilcoxon_fgsea)=c("geneset_id","P.Value","ES")
          result_wilcoxon_fgsea$geneset_id=gsub("\\.", "'", result_wilcoxon_fgsea$geneset_id)
          result_wilcoxon_fgsea=na.omit(result_wilcoxon_fgsea)
          result_wilcoxon_fgsea$two=TRUE
          result_wilcoxon_fgsea$invert_A=ifelse(result_wilcoxon_fgsea$ES<0,TRUE,FALSE)
          result_wilcoxon_fgsea$invert_I=ifelse(result_wilcoxon_fgsea$ES>0,TRUE,FALSE)
          result_wilcoxon_fgsea$P.Value.A=two2one(result_wilcoxon_fgsea$P.Value, two = as.logical(result_wilcoxon_fgsea$two), invert = as.logical(result_wilcoxon_fgsea$invert_A))
          result_wilcoxon_fgsea$P.Value.I=two2one(result_wilcoxon_fgsea$P.Value, two = as.logical(result_wilcoxon_fgsea$two), invert = as.logical(result_wilcoxon_fgsea$invert_I))
          result_wilcoxon_fgsea_A=result_wilcoxon_fgsea[,c("geneset_id","P.Value.A")]
          result_wilcoxon_fgsea_I=result_wilcoxon_fgsea[,c("geneset_id","P.Value.I")]
          colnames(result_wilcoxon_fgsea_A)[2]=paste0(vector[i],".wilcoxon_fgsea")
          colnames(result_wilcoxon_fgsea_I)[2]=paste0(vector[i],".wilcoxon_fgsea")
          rm(result_wilcoxon_fgsea)
        }
        if("wilcoxon_ora" %in% enrichment_method){
          print("wilcoxon_ora will start.")
          result_UP=all_diff[all_diff$logFC>=0,]
          result_DN=all_diff[all_diff$logFC<0,]
          result_UP=result_UP[order(result_UP$P.Value, -result_UP$logFC), ] 
          result_DN=result_DN[order(result_DN$P.Value, result_DN$logFC), ]
          upregulated_gene_sets=head(result_UP,min(1000, nrow(result_UP)))
          downregulated_gene_sets=head(result_DN,min(1000, nrow(result_DN)))
          gene_H=upregulated_gene_sets$geneset_id
          gene_L=downregulated_gene_sets$geneset_id
          canonical_pathways=do.call(rbind, lapply(geneSets, function(gs) {  
            data.frame(term = gs@setName, gene = gs@geneIds, stringsAsFactors = FALSE)  
          }))
          canonical_pathways=canonical_pathways[canonical_pathways$gene != "", ]
          enrich_H=enricher(gene_H, TERM2GENE = canonical_pathways,minGSSize=min.sz,maxGSSize=max.sz,
                            pvalueCutoff = 1,qvalueCutoff = 1)
          if (is.null(enrich_H)) {
            result_wilcoxon_ora_A=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_wilcoxon_ora_A=enrich_H@result
            result_wilcoxon_ora_A=result_wilcoxon_ora_A[,c("Description","pvalue")]
          }
          colnames(result_wilcoxon_ora_A)=c("geneset_id","P.Value_H")
          result_wilcoxon_ora_A$geneset_id=gsub("\\.", "'", result_wilcoxon_ora_A$geneset_id)
          enrich_L=enricher(gene_L, TERM2GENE = canonical_pathways,minGSSize=min.sz,maxGSSize=max.sz,
                            pvalueCutoff = 1,qvalueCutoff = 1)
          if (is.null(enrich_L)) {
            result_wilcoxon_ora_I=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_wilcoxon_ora_I=enrich_L@result
            result_wilcoxon_ora_I=result_wilcoxon_ora_I[,c("Description","pvalue")]
          }
          colnames(result_wilcoxon_ora_I)=c("geneset_id","P.Value_L")
          result_wilcoxon_ora_I$geneset_id=gsub("\\.", "'", result_wilcoxon_ora_I$geneset_id)
          colnames(result_wilcoxon_ora_A)[2]=paste0(vector[i],".wilcoxon_ora")
          colnames(result_wilcoxon_ora_I)[2]=paste0(vector[i],".wilcoxon_ora")
        }
      }
      if (any(c("permutation_fgsea","permutation_ora") %in% enrichment_method)){
        all_diff=Get_permutation_test(geneset=exp,group=group)
        if("permutation_fgsea" %in% enrichment_method){
          print("permutation_fgsea will start.")
          alldiff=all_diff[,c("geneset_id","P.Value","logFC")]
          alldiff$id=-log10(alldiff$P.Value)*alldiff$logFC
          alldiff=alldiff[order(alldiff$id,decreasing = T),]
          id=alldiff$id
          names(id)=alldiff$geneset_id
          canonical_pathways=do.call(rbind, lapply(geneSets, function(gs) {  
            data.frame(term = gs@setName, gene = gs@geneIds, stringsAsFactors = FALSE)  
          }))
          canonical_pathways=canonical_pathways[canonical_pathways$gene != "", ]
          canonical_pathways$term=as.character(canonical_pathways$term)
          canonical_pathways=canonical_pathways %>% split(x = .$gene, f = .$term)
          fgseaRes=fgsea(pathways=canonical_pathways,
                         stats=id[is.finite(id)],
                         minSize=min.sz,
                         maxSize=max.sz,
                         nperm=1000)
          result=data.frame(fgseaRes)
          result_permutation_fgsea=result[,c("pathway","pval","ES")]
          colnames(result_permutation_fgsea)=c("geneset_id","P.Value","ES")
          result_permutation_fgsea$geneset_id=gsub("\\.", "'", result_permutation_fgsea$geneset_id)
          result_permutation_fgsea=na.omit(result_permutation_fgsea)
          result_permutation_fgsea$two=TRUE
          result_permutation_fgsea$invert_A=ifelse(result_permutation_fgsea$ES<0,TRUE,FALSE)
          result_permutation_fgsea$invert_I=ifelse(result_permutation_fgsea$ES>0,TRUE,FALSE)
          result_permutation_fgsea$P.Value.A=two2one(result_permutation_fgsea$P.Value, two = as.logical(result_permutation_fgsea$two), invert = as.logical(result_permutation_fgsea$invert_A))
          result_permutation_fgsea$P.Value.I=two2one(result_permutation_fgsea$P.Value, two = as.logical(result_permutation_fgsea$two), invert = as.logical(result_permutation_fgsea$invert_I))
          result_permutation_fgsea_A=result_permutation_fgsea[,c("geneset_id","P.Value.A")]
          result_permutation_fgsea_I=result_permutation_fgsea[,c("geneset_id","P.Value.I")]
          colnames(result_permutation_fgsea_A)[2]=paste0(vector[i],".permutation_fgsea")
          colnames(result_permutation_fgsea_I)[2]=paste0(vector[i],".permutation_fgsea")
          rm(result_permutation_fgsea)
        }
        if("permutation_ora" %in% enrichment_method){
          print("permutation_ora will start.")
          result_UP=all_diff[all_diff$logFC>=0,]
          result_DN=all_diff[all_diff$logFC<0,]
          result_UP=result_UP[order(result_UP$P.Value, -result_UP$logFC), ] 
          result_DN=result_DN[order(result_DN$P.Value, result_DN$logFC), ]
          upregulated_gene_sets=head(result_UP,min(1000, nrow(result_UP)))
          downregulated_gene_sets=head(result_DN,min(1000, nrow(result_DN)))
          gene_H=upregulated_gene_sets$geneset_id
          gene_L=downregulated_gene_sets$geneset_id
          canonical_pathways=do.call(rbind, lapply(geneSets, function(gs) {  
            data.frame(term = gs@setName, gene = gs@geneIds, stringsAsFactors = FALSE)  
          }))
          canonical_pathways=canonical_pathways[canonical_pathways$gene != "", ]
          enrich_H=enricher(gene_H, TERM2GENE = canonical_pathways,minGSSize=min.sz,maxGSSize=max.sz,
                            pvalueCutoff = 1,qvalueCutoff = 1)
          if (is.null(enrich_H)) {
            result_permutation_ora_A=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_permutation_ora_A=enrich_H@result
            result_permutation_ora_A=result_permutation_ora_A[,c("Description","pvalue")]
          }
          colnames(result_permutation_ora_A)=c("geneset_id","P.Value_H")
          result_permutation_ora_A$geneset_id=gsub("\\.", "'", result_permutation_ora_A$geneset_id)
          enrich_L=enricher(gene_L, TERM2GENE = canonical_pathways,minGSSize=min.sz,maxGSSize=max.sz,
                            pvalueCutoff = 1,qvalueCutoff = 1)
          if (is.null(enrich_L)) {
            result_permutation_ora_I=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_permutation_ora_I=enrich_L@result
            result_permutation_ora_I=result_permutation_ora_I[,c("Description","pvalue")]
          }
          colnames(result_permutation_ora_I)=c("geneset_id","P.Value_L")
          result_permutation_ora_I$geneset_id=gsub("\\.", "'", result_permutation_ora_I$geneset_id)
          colnames(result_permutation_ora_A)[2]=paste0(vector[i],".permutation_ora")
          colnames(result_permutation_ora_I)[2]=paste0(vector[i],".permutation_ora")
        }
      }
      if (any(c("kruskal_fgsea","kruskal_ora") %in% enrichment_method)){
        all_diff=Get_kruskal_test(geneset=exp,group=group)
        if("kruskal_fgsea" %in% enrichment_method){
          print("kruskal_fgsea will start.")
          alldiff=all_diff[,c("geneset_id","P.Value","logFC")]
          alldiff$id=-log10(alldiff$P.Value)*alldiff$logFC
          alldiff=alldiff[order(alldiff$id,decreasing = T),]
          id=alldiff$id
          names(id)=alldiff$geneset_id
          canonical_pathways=do.call(rbind, lapply(geneSets, function(gs) {  
            data.frame(term = gs@setName, gene = gs@geneIds, stringsAsFactors = FALSE)  
          }))
          canonical_pathways=canonical_pathways[canonical_pathways$gene != "", ]
          canonical_pathways$term=as.character(canonical_pathways$term)
          canonical_pathways=canonical_pathways %>% split(x = .$gene, f = .$term)
          fgseaRes=fgsea(pathways=canonical_pathways,
                         stats=id[is.finite(id)],
                         minSize=min.sz,
                         maxSize=max.sz,
                         nperm=1000)
          result=data.frame(fgseaRes)
          result_kruskal_fgsea=result[,c("pathway","pval","ES")]
          colnames(result_kruskal_fgsea)=c("geneset_id","P.Value","ES")
          result_kruskal_fgsea$geneset_id=gsub("\\.", "'", result_kruskal_fgsea$geneset_id)
          result_kruskal_fgsea=na.omit(result_kruskal_fgsea)
          result_kruskal_fgsea$two=TRUE
          result_kruskal_fgsea$invert_A=ifelse(result_kruskal_fgsea$ES<0,TRUE,FALSE)
          result_kruskal_fgsea$invert_I=ifelse(result_kruskal_fgsea$ES>0,TRUE,FALSE)
          result_kruskal_fgsea$P.Value.A=two2one(result_kruskal_fgsea$P.Value, two = as.logical(result_kruskal_fgsea$two), invert = as.logical(result_kruskal_fgsea$invert_A))
          result_kruskal_fgsea$P.Value.I=two2one(result_kruskal_fgsea$P.Value, two = as.logical(result_kruskal_fgsea$two), invert = as.logical(result_kruskal_fgsea$invert_I))
          result_kruskal_fgsea_A=result_kruskal_fgsea[,c("geneset_id","P.Value.A")]
          result_kruskal_fgsea_I=result_kruskal_fgsea[,c("geneset_id","P.Value.I")]
          colnames(result_kruskal_fgsea_A)[2]=paste0(vector[i],".kruskal_fgsea")
          colnames(result_kruskal_fgsea_I)[2]=paste0(vector[i],".kruskal_fgsea")
          rm(result_kruskal_fgsea)
        }
        if("kruskal_ora" %in% enrichment_method){
          print("kruskal_ora will start.")
          result_UP=all_diff[all_diff$logFC>=0,]
          result_DN=all_diff[all_diff$logFC<0,]
          result_UP=result_UP[order(result_UP$P.Value, -result_UP$logFC), ] 
          result_DN=result_DN[order(result_DN$P.Value, result_DN$logFC), ]
          upregulated_gene_sets=head(result_UP,min(1000, nrow(result_UP)))
          downregulated_gene_sets=head(result_DN,min(1000, nrow(result_DN)))
          gene_H=upregulated_gene_sets$geneset_id
          gene_L=downregulated_gene_sets$geneset_id
          canonical_pathways=do.call(rbind, lapply(geneSets, function(gs) {  
            data.frame(term = gs@setName, gene = gs@geneIds, stringsAsFactors = FALSE)  
          }))
          canonical_pathways=canonical_pathways[canonical_pathways$gene != "", ]
          enrich_H=enricher(gene_H, TERM2GENE = canonical_pathways,minGSSize=min.sz,maxGSSize=max.sz,
                            pvalueCutoff = 1,qvalueCutoff = 1)
          if (is.null(enrich_H)) {
            result_kruskal_ora_A=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_kruskal_ora_A=enrich_H@result
            result_kruskal_ora_A=result_kruskal_ora_A[,c("Description","pvalue")]
          }
          colnames(result_kruskal_ora_A)=c("geneset_id","P.Value_H")
          result_kruskal_ora_A$geneset_id=gsub("\\.", "'", result_kruskal_ora_A$geneset_id)
          enrich_L=enricher(gene_L, TERM2GENE = canonical_pathways,minGSSize=min.sz,maxGSSize=max.sz,
                            pvalueCutoff = 1,qvalueCutoff = 1)
          if (is.null(enrich_L)) {
            result_kruskal_ora_I=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_kruskal_ora_I=enrich_L@result
            result_kruskal_ora_I=result_kruskal_ora_I[,c("Description","pvalue")]
          }
          colnames(result_kruskal_ora_I)=c("geneset_id","P.Value_L")
          result_kruskal_ora_I$geneset_id=gsub("\\.", "'", result_kruskal_ora_I$geneset_id)
          colnames(result_kruskal_ora_A)[2]=paste0(vector[i],".kruskal_ora")
          colnames(result_kruskal_ora_I)[2]=paste0(vector[i],".kruskal_ora")
        }
      }
    } else {
      print("Not using fgsea and/or ora")
    }
    ##########################################consensus
    methods_to_check=c("consensus_t", "consensus_limma", "consensus_anova", "consensus_wilcoxon", "consensus_permutation", "consensus_kruskal")
    if (any(methods_to_check %in% enrichment_method)){
      geneSets_list=lapply(geneSets, function(x){ x@geneIds })
      names(geneSets_list)=names(geneSets)
      network=data.frame()
      for(exo in 1:length(geneSets_list)){
        geneSets_list_exo=as.data.frame(geneSets_list[[exo]])
        colnames(geneSets_list_exo)="target"
        geneSets_list_exo$source=names(geneSets_list)[exo]
        geneSets_list_exo=geneSets_list_exo[,c(2,1)]
        geneSets_list_exo$mor=1
        network=rbind(network,geneSets_list_exo)
      }
      network=network[network$target!="",]
      network=network[!grepl("^[0-9]", network$target), ]
      network=unique(network)
      geneset=decouple(mat=exp,network=network,minsize = min.sz)
      geneset=run_consensus(geneset)
      geneset=geneset[,c("condition","source","score")]
      geneset=reshape2::dcast(geneset, source ~ condition, value.var = "score")
      rownames(geneset)=geneset[,1]
      geneset=geneset[,-1]
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("consensus_t" %in% enrichment_method){
        print("consensus_t will start.")
        result=Get_T_test(geneset=geneset,group=group)
        result_consensus_t=result[,c("geneset_id","P.Value","logFC")]
        result_consensus_t$two=TRUE
        result_consensus_t$invert_A=ifelse(result_consensus_t$logFC<=0,TRUE,FALSE)
        result_consensus_t$invert_I=ifelse(result_consensus_t$logFC>0,TRUE,FALSE)
        result_consensus_t$P.Value.A=two2one(result_consensus_t$P.Value, two = as.logical(result_consensus_t$two), invert = as.logical(result_consensus_t$invert_A))
        result_consensus_t$P.Value.I=two2one(result_consensus_t$P.Value, two = as.logical(result_consensus_t$two), invert = as.logical(result_consensus_t$invert_I))
        result_consensus_t_A=result_consensus_t[,c("geneset_id","P.Value.A")]
        result_consensus_t_I=result_consensus_t[,c("geneset_id","P.Value.I")]
        colnames(result_consensus_t_A)[2]=paste0(vector[i],".consensus_t")
        colnames(result_consensus_t_I)[2]=paste0(vector[i],".consensus_t")
        rm(result_consensus_t)
      }
      if("consensus_wilcoxon" %in% enrichment_method){
        print("consensus_wilcoxon will start.")
        result=Get_wilcox_test(geneset=geneset,group=group)
        result_consensus_wilcoxon=result[,c("geneset_id","P.Value","logFC")]
        result_consensus_wilcoxon$two=TRUE
        result_consensus_wilcoxon$invert_A=ifelse(result_consensus_wilcoxon$logFC<=0,TRUE,FALSE)
        result_consensus_wilcoxon$invert_I=ifelse(result_consensus_wilcoxon$logFC>0,TRUE,FALSE)
        result_consensus_wilcoxon$P.Value.A=two2one(result_consensus_wilcoxon$P.Value, two = as.logical(result_consensus_wilcoxon$two), invert = as.logical(result_consensus_wilcoxon$invert_A))
        result_consensus_wilcoxon$P.Value.I=two2one(result_consensus_wilcoxon$P.Value, two = as.logical(result_consensus_wilcoxon$two), invert = as.logical(result_consensus_wilcoxon$invert_I))
        result_consensus_wilcoxon_A=result_consensus_wilcoxon[,c("geneset_id","P.Value.A")]
        result_consensus_wilcoxon_I=result_consensus_wilcoxon[,c("geneset_id","P.Value.I")]
        colnames(result_consensus_wilcoxon_A)[2]=paste0(vector[i],".consensus_wilcoxon")
        colnames(result_consensus_wilcoxon_I)[2]=paste0(vector[i],".consensus_wilcoxon")
        rm(result_consensus_wilcoxon)
      }
      if("consensus_kruskal" %in% enrichment_method){
        print("consensus_kruskal will start.")
        result=Get_kruskal_test(geneset=geneset,group=group)
        result_consensus_kruskal=result[,c("geneset_id","P.Value","logFC")]
        result_consensus_kruskal$two=TRUE
        result_consensus_kruskal$invert_A=ifelse(result_consensus_kruskal$logFC<=0,TRUE,FALSE)
        result_consensus_kruskal$invert_I=ifelse(result_consensus_kruskal$logFC>0,TRUE,FALSE)
        result_consensus_kruskal$P.Value.A=two2one(result_consensus_kruskal$P.Value, two = as.logical(result_consensus_kruskal$two), invert = as.logical(result_consensus_kruskal$invert_A))
        result_consensus_kruskal$P.Value.I=two2one(result_consensus_kruskal$P.Value, two = as.logical(result_consensus_kruskal$two), invert = as.logical(result_consensus_kruskal$invert_I))
        result_consensus_kruskal_A=result_consensus_kruskal[,c("geneset_id","P.Value.A")]
        result_consensus_kruskal_I=result_consensus_kruskal[,c("geneset_id","P.Value.I")]
        colnames(result_consensus_kruskal_A)[2]=paste0(vector[i],".consensus_kruskal")
        colnames(result_consensus_kruskal_I)[2]=paste0(vector[i],".consensus_kruskal")
        rm(result_consensus_kruskal)
      }
      if("consensus_permutation" %in% enrichment_method){
        print("consensus_permutation will start.")
        result=Get_permutation_test(geneset=geneset,group=group)
        result_consensus_permutation=result[,c("geneset_id","P.Value","logFC")]
        result_consensus_permutation$two=TRUE
        result_consensus_permutation$invert_A=ifelse(result_consensus_permutation$logFC<=0,TRUE,FALSE)
        result_consensus_permutation$invert_I=ifelse(result_consensus_permutation$logFC>0,TRUE,FALSE)
        result_consensus_permutation=na.omit(result_consensus_permutation)
        result_consensus_permutation$P.Value.A=two2one(result_consensus_permutation$P.Value, two = as.logical(result_consensus_permutation$two), invert = as.logical(result_consensus_permutation$invert_A))
        result_consensus_permutation$P.Value.I=two2one(result_consensus_permutation$P.Value, two = as.logical(result_consensus_permutation$two), invert = as.logical(result_consensus_permutation$invert_I))
        result_consensus_permutation_A=result_consensus_permutation[,c("geneset_id","P.Value.A")]
        result_consensus_permutation_I=result_consensus_permutation[,c("geneset_id","P.Value.I")]
        colnames(result_consensus_permutation_A)[2]=paste0(vector[i],".consensus_permutation")
        colnames(result_consensus_permutation_I)[2]=paste0(vector[i],".consensus_permutation")
        rm(result_consensus_permutation)
      }
      if("consensus_limma" %in% enrichment_method){
        print("consensus_limma will start.")
        result=Get_limma(geneset=geneset,group=group)
        result_consensus_limma=result[,c("geneset_id","P.Value","logFC")]
        result_consensus_limma$two=TRUE
        result_consensus_limma$invert_A=ifelse(result_consensus_limma$logFC<=0,TRUE,FALSE)
        result_consensus_limma$invert_I=ifelse(result_consensus_limma$logFC>0,TRUE,FALSE)
        result_consensus_limma$P.Value.A=two2one(result_consensus_limma$P.Value, two = as.logical(result_consensus_limma$two), invert = as.logical(result_consensus_limma$invert_A))
        result_consensus_limma$P.Value.I=two2one(result_consensus_limma$P.Value, two = as.logical(result_consensus_limma$two), invert = as.logical(result_consensus_limma$invert_I))
        result_consensus_limma_A=result_consensus_limma[,c("geneset_id","P.Value.A")]
        result_consensus_limma_I=result_consensus_limma[,c("geneset_id","P.Value.I")]
        colnames(result_consensus_limma_A)[2]=paste0(vector[i],".consensus_limma")
        colnames(result_consensus_limma_I)[2]=paste0(vector[i],".consensus_limma")
        rm(result_consensus_limma)
      }
      if("consensus_anova" %in% enrichment_method){
        print("consensus_anova will start.")
        result=Get_anova(geneset=geneset,group=group)
        result_consensus_anova=result[,c("geneset_id","P.Value","logFC")]
        result_consensus_anova$two=TRUE
        result_consensus_anova$invert_A=ifelse(result_consensus_anova$logFC<=0,TRUE,FALSE)
        result_consensus_anova$invert_I=ifelse(result_consensus_anova$logFC>0,TRUE,FALSE)
        result_consensus_anova$P.Value.A=two2one(result_consensus_anova$P.Value, two = as.logical(result_consensus_anova$two), invert = as.logical(result_consensus_anova$invert_A))
        result_consensus_anova$P.Value.I=two2one(result_consensus_anova$P.Value, two = as.logical(result_consensus_anova$two), invert = as.logical(result_consensus_anova$invert_I))
        result_consensus_anova_A=result_consensus_anova[,c("geneset_id","P.Value.A")]
        result_consensus_anova_I=result_consensus_anova[,c("geneset_id","P.Value.I")]
        colnames(result_consensus_anova_A)[2]=paste0(vector[i],".consensus_anova")
        colnames(result_consensus_anova_I)[2]=paste0(vector[i],".consensus_anova")
        rm(result_consensus_anova)
      }
    } else {
      print("Not using consensus.")
    }
    ##########################################mdt
    methods_to_check=c("mdt_t", "mdt_limma", "mdt_anova", "mdt_wilcoxon", "mdt_permutation", "mdt_kruskal")
    if (any(methods_to_check %in% enrichment_method)){
      geneSets_list=lapply(geneSets, function(x){ x@geneIds })
      names(geneSets_list)=names(geneSets)
      network=data.frame()
      for(exo in 1:length(geneSets_list)){
        geneSets_list_exo=as.data.frame(geneSets_list[[exo]])
        colnames(geneSets_list_exo)="target"
        geneSets_list_exo$source=names(geneSets_list)[exo]
        geneSets_list_exo=geneSets_list_exo[,c(2,1)]
        geneSets_list_exo$mor=1
        network=rbind(network,geneSets_list_exo)
      }
      network=network[network$target!="",]
      geneset=run_mdt(mat=exp,network=network,minsize = min.sz)
      geneset=geneset[,c("condition","source","score")]
      geneset=reshape2::dcast(geneset, source ~ condition, value.var = "score")
      rownames(geneset)=geneset[,1]
      geneset=geneset[,-1]
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("mdt_t" %in% enrichment_method){
        print("mdt_t will start.")
        result=Get_T_test(geneset=geneset,group=group)
        result_mdt_t=result[,c("geneset_id","P.Value","logFC")]
        result_mdt_t$two=TRUE
        result_mdt_t$invert_A=ifelse(result_mdt_t$logFC<=0,TRUE,FALSE)
        result_mdt_t$invert_I=ifelse(result_mdt_t$logFC>0,TRUE,FALSE)
        result_mdt_t$P.Value.A=two2one(result_mdt_t$P.Value, two = as.logical(result_mdt_t$two), invert = as.logical(result_mdt_t$invert_A))
        result_mdt_t$P.Value.I=two2one(result_mdt_t$P.Value, two = as.logical(result_mdt_t$two), invert = as.logical(result_mdt_t$invert_I))
        result_mdt_t_A=result_mdt_t[,c("geneset_id","P.Value.A")]
        result_mdt_t_I=result_mdt_t[,c("geneset_id","P.Value.I")]
        colnames(result_mdt_t_A)[2]=paste0(vector[i],".mdt_t")
        colnames(result_mdt_t_I)[2]=paste0(vector[i],".mdt_t")
        rm(result_mdt_t)
      }
      if("mdt_wilcoxon" %in% enrichment_method){
        print("mdt_wilcoxon will start.")
        result=Get_wilcox_test(geneset=geneset,group=group)
        result_mdt_wilcoxon=result[,c("geneset_id","P.Value","logFC")]
        result_mdt_wilcoxon$two=TRUE
        result_mdt_wilcoxon$invert_A=ifelse(result_mdt_wilcoxon$logFC<=0,TRUE,FALSE)
        result_mdt_wilcoxon$invert_I=ifelse(result_mdt_wilcoxon$logFC>0,TRUE,FALSE)
        result_mdt_wilcoxon$P.Value.A=two2one(result_mdt_wilcoxon$P.Value, two = as.logical(result_mdt_wilcoxon$two), invert = as.logical(result_mdt_wilcoxon$invert_A))
        result_mdt_wilcoxon$P.Value.I=two2one(result_mdt_wilcoxon$P.Value, two = as.logical(result_mdt_wilcoxon$two), invert = as.logical(result_mdt_wilcoxon$invert_I))
        result_mdt_wilcoxon_A=result_mdt_wilcoxon[,c("geneset_id","P.Value.A")]
        result_mdt_wilcoxon_I=result_mdt_wilcoxon[,c("geneset_id","P.Value.I")]
        colnames(result_mdt_wilcoxon_A)[2]=paste0(vector[i],".mdt_wilcoxon")
        colnames(result_mdt_wilcoxon_I)[2]=paste0(vector[i],".mdt_wilcoxon")
        rm(result_mdt_wilcoxon)
      }
      if("mdt_kruskal" %in% enrichment_method){
        print("mdt_kruskal will start.")
        result=Get_kruskal_test(geneset=geneset,group=group)
        result_mdt_kruskal=result[,c("geneset_id","P.Value","logFC")]
        result_mdt_kruskal$two=TRUE
        result_mdt_kruskal$invert_A=ifelse(result_mdt_kruskal$logFC<=0,TRUE,FALSE)
        result_mdt_kruskal$invert_I=ifelse(result_mdt_kruskal$logFC>0,TRUE,FALSE)
        result_mdt_kruskal$P.Value.A=two2one(result_mdt_kruskal$P.Value, two = as.logical(result_mdt_kruskal$two), invert = as.logical(result_mdt_kruskal$invert_A))
        result_mdt_kruskal$P.Value.I=two2one(result_mdt_kruskal$P.Value, two = as.logical(result_mdt_kruskal$two), invert = as.logical(result_mdt_kruskal$invert_I))
        result_mdt_kruskal_A=result_mdt_kruskal[,c("geneset_id","P.Value.A")]
        result_mdt_kruskal_I=result_mdt_kruskal[,c("geneset_id","P.Value.I")]
        colnames(result_mdt_kruskal_A)[2]=paste0(vector[i],".mdt_kruskal")
        colnames(result_mdt_kruskal_I)[2]=paste0(vector[i],".mdt_kruskal")
        rm(result_mdt_kruskal)
      }
      if("mdt_permutation" %in% enrichment_method){
        print("mdt_permutation will start.")
        result=Get_permutation_test(geneset=geneset,group=group)
        result_mdt_permutation=result[,c("geneset_id","P.Value","logFC")]
        result_mdt_permutation$two=TRUE
        result_mdt_permutation$invert_A=ifelse(result_mdt_permutation$logFC<=0,TRUE,FALSE)
        result_mdt_permutation$invert_I=ifelse(result_mdt_permutation$logFC>0,TRUE,FALSE)
        result_mdt_permutation=na.omit(result_mdt_permutation)
        result_mdt_permutation$P.Value.A=two2one(result_mdt_permutation$P.Value, two = as.logical(result_mdt_permutation$two), invert = as.logical(result_mdt_permutation$invert_A))
        result_mdt_permutation$P.Value.I=two2one(result_mdt_permutation$P.Value, two = as.logical(result_mdt_permutation$two), invert = as.logical(result_mdt_permutation$invert_I))
        result_mdt_permutation_A=result_mdt_permutation[,c("geneset_id","P.Value.A")]
        result_mdt_permutation_I=result_mdt_permutation[,c("geneset_id","P.Value.I")]
        colnames(result_mdt_permutation_A)[2]=paste0(vector[i],".mdt_permutation")
        colnames(result_mdt_permutation_I)[2]=paste0(vector[i],".mdt_permutation")
        rm(result_mdt_permutation)
      }
      if("mdt_limma" %in% enrichment_method){
        print("mdt_limma will start.")
        result=Get_limma(geneset=geneset,group=group)
        result_mdt_limma=result[,c("geneset_id","P.Value","logFC")]
        result_mdt_limma$two=TRUE
        result_mdt_limma$invert_A=ifelse(result_mdt_limma$logFC<=0,TRUE,FALSE)
        result_mdt_limma$invert_I=ifelse(result_mdt_limma$logFC>0,TRUE,FALSE)
        result_mdt_limma$P.Value.A=two2one(result_mdt_limma$P.Value, two = as.logical(result_mdt_limma$two), invert = as.logical(result_mdt_limma$invert_A))
        result_mdt_limma$P.Value.I=two2one(result_mdt_limma$P.Value, two = as.logical(result_mdt_limma$two), invert = as.logical(result_mdt_limma$invert_I))
        result_mdt_limma_A=result_mdt_limma[,c("geneset_id","P.Value.A")]
        result_mdt_limma_I=result_mdt_limma[,c("geneset_id","P.Value.I")]
        colnames(result_mdt_limma_A)[2]=paste0(vector[i],".mdt_limma")
        colnames(result_mdt_limma_I)[2]=paste0(vector[i],".mdt_limma")
        rm(result_mdt_limma)
      }
      if("mdt_anova" %in% enrichment_method){
        print("mdt_anova will start.")
        result=Get_anova(geneset=geneset,group=group)
        result_mdt_anova=result[,c("geneset_id","P.Value","logFC")]
        result_mdt_anova$two=TRUE
        result_mdt_anova$invert_A=ifelse(result_mdt_anova$logFC<=0,TRUE,FALSE)
        result_mdt_anova$invert_I=ifelse(result_mdt_anova$logFC>0,TRUE,FALSE)
        result_mdt_anova$P.Value.A=two2one(result_mdt_anova$P.Value, two = as.logical(result_mdt_anova$two), invert = as.logical(result_mdt_anova$invert_A))
        result_mdt_anova$P.Value.I=two2one(result_mdt_anova$P.Value, two = as.logical(result_mdt_anova$two), invert = as.logical(result_mdt_anova$invert_I))
        result_mdt_anova_A=result_mdt_anova[,c("geneset_id","P.Value.A")]
        result_mdt_anova_I=result_mdt_anova[,c("geneset_id","P.Value.I")]
        colnames(result_mdt_anova_A)[2]=paste0(vector[i],".mdt_anova")
        colnames(result_mdt_anova_I)[2]=paste0(vector[i],".mdt_anova")
        rm(result_mdt_anova)
      }
    } else {
      print("Not using mdt.")
    }
    ##########################################mlm
    methods_to_check=c("mlm_t", "mlm_limma", "mlm_anova", "mlm_wilcoxon", "mlm_permutation", "mlm_kruskal")
    if (any(methods_to_check %in% enrichment_method)){
      geneSets_list=lapply(geneSets, function(x){ x@geneIds })
      names(geneSets_list)=names(geneSets)
      network=data.frame()
      for(exo in 1:length(geneSets_list)){
        geneSets_list_exo=as.data.frame(geneSets_list[[exo]])
        colnames(geneSets_list_exo)="target"
        geneSets_list_exo$source=names(geneSets_list)[exo]
        geneSets_list_exo=geneSets_list_exo[,c(2,1)]
        geneSets_list_exo$mor=1
        network=rbind(network,geneSets_list_exo)
      }
      network=network[network$target!="",]
      geneset=run_mlm(mat=exp,network=network,minsize = min.sz)
      geneset=geneset[,c("condition","source","score")]
      geneset=reshape2::dcast(geneset, source ~ condition, value.var = "score")
      rownames(geneset)=geneset[,1]
      geneset=geneset[,-1]
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("mlm_t" %in% enrichment_method){
        print("mlm_t will start.")
        result=Get_T_test(geneset=geneset,group=group)
        result_mlm_t=result[,c("geneset_id","P.Value","logFC")]
        result_mlm_t$two=TRUE
        result_mlm_t$invert_A=ifelse(result_mlm_t$logFC<=0,TRUE,FALSE)
        result_mlm_t$invert_I=ifelse(result_mlm_t$logFC>0,TRUE,FALSE)
        result_mlm_t$P.Value.A=two2one(result_mlm_t$P.Value, two = as.logical(result_mlm_t$two), invert = as.logical(result_mlm_t$invert_A))
        result_mlm_t$P.Value.I=two2one(result_mlm_t$P.Value, two = as.logical(result_mlm_t$two), invert = as.logical(result_mlm_t$invert_I))
        result_mlm_t_A=result_mlm_t[,c("geneset_id","P.Value.A")]
        result_mlm_t_I=result_mlm_t[,c("geneset_id","P.Value.I")]
        colnames(result_mlm_t_A)[2]=paste0(vector[i],".mlm_t")
        colnames(result_mlm_t_I)[2]=paste0(vector[i],".mlm_t")
        rm(result_mlm_t)
      }
      if("mlm_wilcoxon" %in% enrichment_method){
        print("mlm_wilcoxon will start.")
        result=Get_wilcox_test(geneset=geneset,group=group)
        result_mlm_wilcoxon=result[,c("geneset_id","P.Value","logFC")]
        result_mlm_wilcoxon$two=TRUE
        result_mlm_wilcoxon$invert_A=ifelse(result_mlm_wilcoxon$logFC<=0,TRUE,FALSE)
        result_mlm_wilcoxon$invert_I=ifelse(result_mlm_wilcoxon$logFC>0,TRUE,FALSE)
        result_mlm_wilcoxon$P.Value.A=two2one(result_mlm_wilcoxon$P.Value, two = as.logical(result_mlm_wilcoxon$two), invert = as.logical(result_mlm_wilcoxon$invert_A))
        result_mlm_wilcoxon$P.Value.I=two2one(result_mlm_wilcoxon$P.Value, two = as.logical(result_mlm_wilcoxon$two), invert = as.logical(result_mlm_wilcoxon$invert_I))
        result_mlm_wilcoxon_A=result_mlm_wilcoxon[,c("geneset_id","P.Value.A")]
        result_mlm_wilcoxon_I=result_mlm_wilcoxon[,c("geneset_id","P.Value.I")]
        colnames(result_mlm_wilcoxon_A)[2]=paste0(vector[i],".mlm_wilcoxon")
        colnames(result_mlm_wilcoxon_I)[2]=paste0(vector[i],".mlm_wilcoxon")
        rm(result_mlm_wilcoxon)
      }
      if("mlm_kruskal" %in% enrichment_method){
        print("mlm_kruskal will start.")
        result=Get_kruskal_test(geneset=geneset,group=group)
        result_mlm_kruskal=result[,c("geneset_id","P.Value","logFC")]
        result_mlm_kruskal$two=TRUE
        result_mlm_kruskal$invert_A=ifelse(result_mlm_kruskal$logFC<=0,TRUE,FALSE)
        result_mlm_kruskal$invert_I=ifelse(result_mlm_kruskal$logFC>0,TRUE,FALSE)
        result_mlm_kruskal$P.Value.A=two2one(result_mlm_kruskal$P.Value, two = as.logical(result_mlm_kruskal$two), invert = as.logical(result_mlm_kruskal$invert_A))
        result_mlm_kruskal$P.Value.I=two2one(result_mlm_kruskal$P.Value, two = as.logical(result_mlm_kruskal$two), invert = as.logical(result_mlm_kruskal$invert_I))
        result_mlm_kruskal_A=result_mlm_kruskal[,c("geneset_id","P.Value.A")]
        result_mlm_kruskal_I=result_mlm_kruskal[,c("geneset_id","P.Value.I")]
        colnames(result_mlm_kruskal_A)[2]=paste0(vector[i],".mlm_kruskal")
        colnames(result_mlm_kruskal_I)[2]=paste0(vector[i],".mlm_kruskal")
        rm(result_mlm_kruskal)
      }
      if("mlm_permutation" %in% enrichment_method){
        print("mlm_permutation will start.")
        result=Get_permutation_test(geneset=geneset,group=group)
        result_mlm_permutation=result[,c("geneset_id","P.Value","logFC")]
        result_mlm_permutation$two=TRUE
        result_mlm_permutation$invert_A=ifelse(result_mlm_permutation$logFC<=0,TRUE,FALSE)
        result_mlm_permutation$invert_I=ifelse(result_mlm_permutation$logFC>0,TRUE,FALSE)
        result_mlm_permutation=na.omit(result_mlm_permutation)
        result_mlm_permutation$P.Value.A=two2one(result_mlm_permutation$P.Value, two = as.logical(result_mlm_permutation$two), invert = as.logical(result_mlm_permutation$invert_A))
        result_mlm_permutation$P.Value.I=two2one(result_mlm_permutation$P.Value, two = as.logical(result_mlm_permutation$two), invert = as.logical(result_mlm_permutation$invert_I))
        result_mlm_permutation_A=result_mlm_permutation[,c("geneset_id","P.Value.A")]
        result_mlm_permutation_I=result_mlm_permutation[,c("geneset_id","P.Value.I")]
        colnames(result_mlm_permutation_A)[2]=paste0(vector[i],".mlm_permutation")
        colnames(result_mlm_permutation_I)[2]=paste0(vector[i],".mlm_permutation")
        rm(result_mlm_permutation)
      }
      if("mlm_limma" %in% enrichment_method){
        print("mlm_limma will start.")
        result=Get_limma(geneset=geneset,group=group)
        result_mlm_limma=result[,c("geneset_id","P.Value","logFC")]
        result_mlm_limma$two=TRUE
        result_mlm_limma$invert_A=ifelse(result_mlm_limma$logFC<=0,TRUE,FALSE)
        result_mlm_limma$invert_I=ifelse(result_mlm_limma$logFC>0,TRUE,FALSE)
        result_mlm_limma$P.Value.A=two2one(result_mlm_limma$P.Value, two = as.logical(result_mlm_limma$two), invert = as.logical(result_mlm_limma$invert_A))
        result_mlm_limma$P.Value.I=two2one(result_mlm_limma$P.Value, two = as.logical(result_mlm_limma$two), invert = as.logical(result_mlm_limma$invert_I))
        result_mlm_limma_A=result_mlm_limma[,c("geneset_id","P.Value.A")]
        result_mlm_limma_I=result_mlm_limma[,c("geneset_id","P.Value.I")]
        colnames(result_mlm_limma_A)[2]=paste0(vector[i],".mlm_limma")
        colnames(result_mlm_limma_I)[2]=paste0(vector[i],".mlm_limma")
        rm(result_mlm_limma)
      }
      if("mlm_anova" %in% enrichment_method){
        print("mlm_anova will start.")
        result=Get_anova(geneset=geneset,group=group)
        result_mlm_anova=result[,c("geneset_id","P.Value","logFC")]
        result_mlm_anova$two=TRUE
        result_mlm_anova$invert_A=ifelse(result_mlm_anova$logFC<=0,TRUE,FALSE)
        result_mlm_anova$invert_I=ifelse(result_mlm_anova$logFC>0,TRUE,FALSE)
        result_mlm_anova$P.Value.A=two2one(result_mlm_anova$P.Value, two = as.logical(result_mlm_anova$two), invert = as.logical(result_mlm_anova$invert_A))
        result_mlm_anova$P.Value.I=two2one(result_mlm_anova$P.Value, two = as.logical(result_mlm_anova$two), invert = as.logical(result_mlm_anova$invert_I))
        result_mlm_anova_A=result_mlm_anova[,c("geneset_id","P.Value.A")]
        result_mlm_anova_I=result_mlm_anova[,c("geneset_id","P.Value.I")]
        colnames(result_mlm_anova_A)[2]=paste0(vector[i],".mlm_anova")
        colnames(result_mlm_anova_I)[2]=paste0(vector[i],".mlm_anova")
        rm(result_mlm_anova)
      }
    } else {
      print("Not using mlm.")
    }
    ##########################################udt
    methods_to_check=c("udt_t", "udt_limma", "udt_anova", "udt_wilcoxon", "udt_permutation", "udt_kruskal")
    if (any(methods_to_check %in% enrichment_method)){
      geneSets_list=lapply(geneSets, function(x){ x@geneIds })
      names(geneSets_list)=names(geneSets)
      network=data.frame()
      for(exo in 1:length(geneSets_list)){
        geneSets_list_exo=as.data.frame(geneSets_list[[exo]])
        colnames(geneSets_list_exo)="target"
        geneSets_list_exo$source=names(geneSets_list)[exo]
        geneSets_list_exo=geneSets_list_exo[,c(2,1)]
        geneSets_list_exo$mor=1
        network=rbind(network,geneSets_list_exo)
      }
      network=network[network$target!="",]
      geneset=run_udt(mat=exp,network=network,minsize = min.sz)
      geneset=geneset[,c("condition","source","score")]
      geneset=reshape2::dcast(geneset, source ~ condition, value.var = "score")
      rownames(geneset)=geneset[,1]
      geneset=geneset[,-1]
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("udt_t" %in% enrichment_method){
        print("udt_t will start.")
        result=Get_T_test(geneset=geneset,group=group)
        result_udt_t=result[,c("geneset_id","P.Value","logFC")]
        result_udt_t$two=TRUE
        result_udt_t$invert_A=ifelse(result_udt_t$logFC<=0,TRUE,FALSE)
        result_udt_t$invert_I=ifelse(result_udt_t$logFC>0,TRUE,FALSE)
        result_udt_t$P.Value.A=two2one(result_udt_t$P.Value, two = as.logical(result_udt_t$two), invert = as.logical(result_udt_t$invert_A))
        result_udt_t$P.Value.I=two2one(result_udt_t$P.Value, two = as.logical(result_udt_t$two), invert = as.logical(result_udt_t$invert_I))
        result_udt_t_A=result_udt_t[,c("geneset_id","P.Value.A")]
        result_udt_t_I=result_udt_t[,c("geneset_id","P.Value.I")]
        colnames(result_udt_t_A)[2]=paste0(vector[i],".udt_t")
        colnames(result_udt_t_I)[2]=paste0(vector[i],".udt_t")
        rm(result_udt_t)
      }
      if("udt_wilcoxon" %in% enrichment_method){
        print("udt_wilcoxon will start.")
        result=Get_wilcox_test(geneset=geneset,group=group)
        result_udt_wilcoxon=result[,c("geneset_id","P.Value","logFC")]
        result_udt_wilcoxon$two=TRUE
        result_udt_wilcoxon$invert_A=ifelse(result_udt_wilcoxon$logFC<=0,TRUE,FALSE)
        result_udt_wilcoxon$invert_I=ifelse(result_udt_wilcoxon$logFC>0,TRUE,FALSE)
        result_udt_wilcoxon$P.Value.A=two2one(result_udt_wilcoxon$P.Value, two = as.logical(result_udt_wilcoxon$two), invert = as.logical(result_udt_wilcoxon$invert_A))
        result_udt_wilcoxon$P.Value.I=two2one(result_udt_wilcoxon$P.Value, two = as.logical(result_udt_wilcoxon$two), invert = as.logical(result_udt_wilcoxon$invert_I))
        result_udt_wilcoxon_A=result_udt_wilcoxon[,c("geneset_id","P.Value.A")]
        result_udt_wilcoxon_I=result_udt_wilcoxon[,c("geneset_id","P.Value.I")]
        colnames(result_udt_wilcoxon_A)[2]=paste0(vector[i],".udt_wilcoxon")
        colnames(result_udt_wilcoxon_I)[2]=paste0(vector[i],".udt_wilcoxon")
        rm(result_udt_wilcoxon)
      }
      if("udt_kruskal" %in% enrichment_method){
        print("udt_kruskal will start.")
        result=Get_kruskal_test(geneset=geneset,group=group)
        result_udt_kruskal=result[,c("geneset_id","P.Value","logFC")]
        result_udt_kruskal$two=TRUE
        result_udt_kruskal$invert_A=ifelse(result_udt_kruskal$logFC<=0,TRUE,FALSE)
        result_udt_kruskal$invert_I=ifelse(result_udt_kruskal$logFC>0,TRUE,FALSE)
        result_udt_kruskal$P.Value.A=two2one(result_udt_kruskal$P.Value, two = as.logical(result_udt_kruskal$two), invert = as.logical(result_udt_kruskal$invert_A))
        result_udt_kruskal$P.Value.I=two2one(result_udt_kruskal$P.Value, two = as.logical(result_udt_kruskal$two), invert = as.logical(result_udt_kruskal$invert_I))
        result_udt_kruskal_A=result_udt_kruskal[,c("geneset_id","P.Value.A")]
        result_udt_kruskal_I=result_udt_kruskal[,c("geneset_id","P.Value.I")]
        colnames(result_udt_kruskal_A)[2]=paste0(vector[i],".udt_kruskal")
        colnames(result_udt_kruskal_I)[2]=paste0(vector[i],".udt_kruskal")
        rm(result_udt_kruskal)
      }
      if("udt_permutation" %in% enrichment_method){
        print("udt_permutation will start.")
        result=Get_permutation_test(geneset=geneset,group=group)
        result_udt_permutation=result[,c("geneset_id","P.Value","logFC")]
        result_udt_permutation$two=TRUE
        result_udt_permutation$invert_A=ifelse(result_udt_permutation$logFC<=0,TRUE,FALSE)
        result_udt_permutation$invert_I=ifelse(result_udt_permutation$logFC>0,TRUE,FALSE)
        result_udt_permutation=na.omit(result_udt_permutation)
        result_udt_permutation$P.Value.A=two2one(result_udt_permutation$P.Value, two = as.logical(result_udt_permutation$two), invert = as.logical(result_udt_permutation$invert_A))
        result_udt_permutation$P.Value.I=two2one(result_udt_permutation$P.Value, two = as.logical(result_udt_permutation$two), invert = as.logical(result_udt_permutation$invert_I))
        result_udt_permutation_A=result_udt_permutation[,c("geneset_id","P.Value.A")]
        result_udt_permutation_I=result_udt_permutation[,c("geneset_id","P.Value.I")]
        colnames(result_udt_permutation_A)[2]=paste0(vector[i],".udt_permutation")
        colnames(result_udt_permutation_I)[2]=paste0(vector[i],".udt_permutation")
        rm(result_udt_permutation)
      }
      if("udt_limma" %in% enrichment_method){
        print("udt_limma will start.")
        result=Get_limma(geneset=geneset,group=group)
        result_udt_limma=result[,c("geneset_id","P.Value","logFC")]
        result_udt_limma$two=TRUE
        result_udt_limma$invert_A=ifelse(result_udt_limma$logFC<=0,TRUE,FALSE)
        result_udt_limma$invert_I=ifelse(result_udt_limma$logFC>0,TRUE,FALSE)
        result_udt_limma$P.Value.A=two2one(result_udt_limma$P.Value, two = as.logical(result_udt_limma$two), invert = as.logical(result_udt_limma$invert_A))
        result_udt_limma$P.Value.I=two2one(result_udt_limma$P.Value, two = as.logical(result_udt_limma$two), invert = as.logical(result_udt_limma$invert_I))
        result_udt_limma_A=result_udt_limma[,c("geneset_id","P.Value.A")]
        result_udt_limma_I=result_udt_limma[,c("geneset_id","P.Value.I")]
        colnames(result_udt_limma_A)[2]=paste0(vector[i],".udt_limma")
        colnames(result_udt_limma_I)[2]=paste0(vector[i],".udt_limma")
        rm(result_udt_limma)
      }
      if("udt_anova" %in% enrichment_method){
        print("udt_anova will start.")
        result=Get_anova(geneset=geneset,group=group)
        result_udt_anova=result[,c("geneset_id","P.Value","logFC")]
        result_udt_anova$two=TRUE
        result_udt_anova$invert_A=ifelse(result_udt_anova$logFC<=0,TRUE,FALSE)
        result_udt_anova$invert_I=ifelse(result_udt_anova$logFC>0,TRUE,FALSE)
        result_udt_anova$P.Value.A=two2one(result_udt_anova$P.Value, two = as.logical(result_udt_anova$two), invert = as.logical(result_udt_anova$invert_A))
        result_udt_anova$P.Value.I=two2one(result_udt_anova$P.Value, two = as.logical(result_udt_anova$two), invert = as.logical(result_udt_anova$invert_I))
        result_udt_anova_A=result_udt_anova[,c("geneset_id","P.Value.A")]
        result_udt_anova_I=result_udt_anova[,c("geneset_id","P.Value.I")]
        colnames(result_udt_anova_A)[2]=paste0(vector[i],".udt_anova")
        colnames(result_udt_anova_I)[2]=paste0(vector[i],".udt_anova")
        rm(result_udt_anova)
      }
    } else {
      print("Not using udt.")
    }
    ##########################################ulm
    methods_to_check=c("ulm_t", "ulm_limma", "ulm_anova", "ulm_wilcoxon", "ulm_permutation", "ulm_kruskal")
    if (any(methods_to_check %in% enrichment_method)){
      geneSets_list=lapply(geneSets, function(x){ x@geneIds })
      names(geneSets_list)=names(geneSets)
      network=data.frame()
      for(exo in 1:length(geneSets_list)){
        geneSets_list_exo=as.data.frame(geneSets_list[[exo]])
        colnames(geneSets_list_exo)="target"
        geneSets_list_exo$source=names(geneSets_list)[exo]
        geneSets_list_exo=geneSets_list_exo[,c(2,1)]
        geneSets_list_exo$mor=1
        network=rbind(network,geneSets_list_exo)
      }
      network=network[network$target!="",]
      geneset=run_ulm(mat=exp,network=network,minsize = min.sz)
      geneset=geneset[,c("condition","source","score")]
      geneset=reshape2::dcast(geneset, source ~ condition, value.var = "score")
      rownames(geneset)=geneset[,1]
      geneset=geneset[,-1]
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("ulm_t" %in% enrichment_method){
        print("ulm_t will start.")
        result=Get_T_test(geneset=geneset,group=group)
        result_ulm_t=result[,c("geneset_id","P.Value","logFC")]
        result_ulm_t$two=TRUE
        result_ulm_t$invert_A=ifelse(result_ulm_t$logFC<=0,TRUE,FALSE)
        result_ulm_t$invert_I=ifelse(result_ulm_t$logFC>0,TRUE,FALSE)
        result_ulm_t$P.Value.A=two2one(result_ulm_t$P.Value, two = as.logical(result_ulm_t$two), invert = as.logical(result_ulm_t$invert_A))
        result_ulm_t$P.Value.I=two2one(result_ulm_t$P.Value, two = as.logical(result_ulm_t$two), invert = as.logical(result_ulm_t$invert_I))
        result_ulm_t_A=result_ulm_t[,c("geneset_id","P.Value.A")]
        result_ulm_t_I=result_ulm_t[,c("geneset_id","P.Value.I")]
        colnames(result_ulm_t_A)[2]=paste0(vector[i],".ulm_t")
        colnames(result_ulm_t_I)[2]=paste0(vector[i],".ulm_t")
        rm(result_ulm_t)
      }
      if("ulm_wilcoxon" %in% enrichment_method){
        print("ulm_wilcoxon will start.")
        result=Get_wilcox_test(geneset=geneset,group=group)
        result_ulm_wilcoxon=result[,c("geneset_id","P.Value","logFC")]
        result_ulm_wilcoxon$two=TRUE
        result_ulm_wilcoxon$invert_A=ifelse(result_ulm_wilcoxon$logFC<=0,TRUE,FALSE)
        result_ulm_wilcoxon$invert_I=ifelse(result_ulm_wilcoxon$logFC>0,TRUE,FALSE)
        result_ulm_wilcoxon$P.Value.A=two2one(result_ulm_wilcoxon$P.Value, two = as.logical(result_ulm_wilcoxon$two), invert = as.logical(result_ulm_wilcoxon$invert_A))
        result_ulm_wilcoxon$P.Value.I=two2one(result_ulm_wilcoxon$P.Value, two = as.logical(result_ulm_wilcoxon$two), invert = as.logical(result_ulm_wilcoxon$invert_I))
        result_ulm_wilcoxon_A=result_ulm_wilcoxon[,c("geneset_id","P.Value.A")]
        result_ulm_wilcoxon_I=result_ulm_wilcoxon[,c("geneset_id","P.Value.I")]
        colnames(result_ulm_wilcoxon_A)[2]=paste0(vector[i],".ulm_wilcoxon")
        colnames(result_ulm_wilcoxon_I)[2]=paste0(vector[i],".ulm_wilcoxon")
        rm(result_ulm_wilcoxon)
      }
      if("ulm_kruskal" %in% enrichment_method){
        print("ulm_kruskal will start.")
        result=Get_kruskal_test(geneset=geneset,group=group)
        result_ulm_kruskal=result[,c("geneset_id","P.Value","logFC")]
        result_ulm_kruskal$two=TRUE
        result_ulm_kruskal$invert_A=ifelse(result_ulm_kruskal$logFC<=0,TRUE,FALSE)
        result_ulm_kruskal$invert_I=ifelse(result_ulm_kruskal$logFC>0,TRUE,FALSE)
        result_ulm_kruskal$P.Value.A=two2one(result_ulm_kruskal$P.Value, two = as.logical(result_ulm_kruskal$two), invert = as.logical(result_ulm_kruskal$invert_A))
        result_ulm_kruskal$P.Value.I=two2one(result_ulm_kruskal$P.Value, two = as.logical(result_ulm_kruskal$two), invert = as.logical(result_ulm_kruskal$invert_I))
        result_ulm_kruskal_A=result_ulm_kruskal[,c("geneset_id","P.Value.A")]
        result_ulm_kruskal_I=result_ulm_kruskal[,c("geneset_id","P.Value.I")]
        colnames(result_ulm_kruskal_A)[2]=paste0(vector[i],".ulm_kruskal")
        colnames(result_ulm_kruskal_I)[2]=paste0(vector[i],".ulm_kruskal")
        rm(result_ulm_kruskal)
      }
      if("ulm_permutation" %in% enrichment_method){
        print("ulm_permutation will start.")
        result=Get_permutation_test(geneset=geneset,group=group)
        result_ulm_permutation=result[,c("geneset_id","P.Value","logFC")]
        result_ulm_permutation$two=TRUE
        result_ulm_permutation$invert_A=ifelse(result_ulm_permutation$logFC<=0,TRUE,FALSE)
        result_ulm_permutation$invert_I=ifelse(result_ulm_permutation$logFC>0,TRUE,FALSE)
        result_ulm_permutation=na.omit(result_ulm_permutation)
        result_ulm_permutation$P.Value.A=two2one(result_ulm_permutation$P.Value, two = as.logical(result_ulm_permutation$two), invert = as.logical(result_ulm_permutation$invert_A))
        result_ulm_permutation$P.Value.I=two2one(result_ulm_permutation$P.Value, two = as.logical(result_ulm_permutation$two), invert = as.logical(result_ulm_permutation$invert_I))
        result_ulm_permutation_A=result_ulm_permutation[,c("geneset_id","P.Value.A")]
        result_ulm_permutation_I=result_ulm_permutation[,c("geneset_id","P.Value.I")]
        colnames(result_ulm_permutation_A)[2]=paste0(vector[i],".ulm_permutation")
        colnames(result_ulm_permutation_I)[2]=paste0(vector[i],".ulm_permutation")
        rm(result_ulm_permutation)
      }
      if("ulm_limma" %in% enrichment_method){
        print("ulm_limma will start.")
        result=Get_limma(geneset=geneset,group=group)
        result_ulm_limma=result[,c("geneset_id","P.Value","logFC")]
        result_ulm_limma$two=TRUE
        result_ulm_limma$invert_A=ifelse(result_ulm_limma$logFC<=0,TRUE,FALSE)
        result_ulm_limma$invert_I=ifelse(result_ulm_limma$logFC>0,TRUE,FALSE)
        result_ulm_limma$P.Value.A=two2one(result_ulm_limma$P.Value, two = as.logical(result_ulm_limma$two), invert = as.logical(result_ulm_limma$invert_A))
        result_ulm_limma$P.Value.I=two2one(result_ulm_limma$P.Value, two = as.logical(result_ulm_limma$two), invert = as.logical(result_ulm_limma$invert_I))
        result_ulm_limma_A=result_ulm_limma[,c("geneset_id","P.Value.A")]
        result_ulm_limma_I=result_ulm_limma[,c("geneset_id","P.Value.I")]
        colnames(result_ulm_limma_A)[2]=paste0(vector[i],".ulm_limma")
        colnames(result_ulm_limma_I)[2]=paste0(vector[i],".ulm_limma")
        rm(result_ulm_limma)
      }
      if("ulm_anova" %in% enrichment_method){
        print("ulm_anova will start.")
        result=Get_anova(geneset=geneset,group=group)
        result_ulm_anova=result[,c("geneset_id","P.Value","logFC")]
        result_ulm_anova$two=TRUE
        result_ulm_anova$invert_A=ifelse(result_ulm_anova$logFC<=0,TRUE,FALSE)
        result_ulm_anova$invert_I=ifelse(result_ulm_anova$logFC>0,TRUE,FALSE)
        result_ulm_anova$P.Value.A=two2one(result_ulm_anova$P.Value, two = as.logical(result_ulm_anova$two), invert = as.logical(result_ulm_anova$invert_A))
        result_ulm_anova$P.Value.I=two2one(result_ulm_anova$P.Value, two = as.logical(result_ulm_anova$two), invert = as.logical(result_ulm_anova$invert_I))
        result_ulm_anova_A=result_ulm_anova[,c("geneset_id","P.Value.A")]
        result_ulm_anova_I=result_ulm_anova[,c("geneset_id","P.Value.I")]
        colnames(result_ulm_anova_A)[2]=paste0(vector[i],".ulm_anova")
        colnames(result_ulm_anova_I)[2]=paste0(vector[i],".ulm_anova")
        rm(result_ulm_anova)
      }
    } else {
      print("Not using ulm.")
    }
    ##########################################viper
    methods_to_check=c("viper_t", "viper_limma", "viper_anova", "viper_wilcoxon", "viper_permutation", "viper_kruskal")
    if (any(methods_to_check %in% enrichment_method)){
      geneSets_list=lapply(geneSets, function(x){ x@geneIds })
      names(geneSets_list)=names(geneSets)
      network=data.frame()
      for(exo in 1:length(geneSets_list)){
        geneSets_list_exo=as.data.frame(geneSets_list[[exo]])
        colnames(geneSets_list_exo)="target"
        geneSets_list_exo$source=names(geneSets_list)[exo]
        geneSets_list_exo=geneSets_list_exo[,c(2,1)]
        geneSets_list_exo$mor=1
        network=rbind(network,geneSets_list_exo)
      }
      network=network[network$target!="",]
      geneset=run_viper(mat=exp,network=network,minsize = min.sz)
      geneset=geneset[,c("condition","source","score")]
      geneset=reshape2::dcast(geneset, source ~ condition, value.var = "score")
      rownames(geneset)=geneset[,1]
      geneset=geneset[,-1]
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("viper_t" %in% enrichment_method){
        print("viper_t will start.")
        result=Get_T_test(geneset=geneset,group=group)
        result_viper_t=result[,c("geneset_id","P.Value","logFC")]
        result_viper_t$two=TRUE
        result_viper_t$invert_A=ifelse(result_viper_t$logFC<=0,TRUE,FALSE)
        result_viper_t$invert_I=ifelse(result_viper_t$logFC>0,TRUE,FALSE)
        result_viper_t$P.Value.A=two2one(result_viper_t$P.Value, two = as.logical(result_viper_t$two), invert = as.logical(result_viper_t$invert_A))
        result_viper_t$P.Value.I=two2one(result_viper_t$P.Value, two = as.logical(result_viper_t$two), invert = as.logical(result_viper_t$invert_I))
        result_viper_t_A=result_viper_t[,c("geneset_id","P.Value.A")]
        result_viper_t_I=result_viper_t[,c("geneset_id","P.Value.I")]
        colnames(result_viper_t_A)[2]=paste0(vector[i],".viper_t")
        colnames(result_viper_t_I)[2]=paste0(vector[i],".viper_t")
        rm(result_viper_t)
      }
      if("viper_wilcoxon" %in% enrichment_method){
        print("viper_wilcoxon will start.")
        result=Get_wilcox_test(geneset=geneset,group=group)
        result_viper_wilcoxon=result[,c("geneset_id","P.Value","logFC")]
        result_viper_wilcoxon$two=TRUE
        result_viper_wilcoxon$invert_A=ifelse(result_viper_wilcoxon$logFC<=0,TRUE,FALSE)
        result_viper_wilcoxon$invert_I=ifelse(result_viper_wilcoxon$logFC>0,TRUE,FALSE)
        result_viper_wilcoxon$P.Value.A=two2one(result_viper_wilcoxon$P.Value, two = as.logical(result_viper_wilcoxon$two), invert = as.logical(result_viper_wilcoxon$invert_A))
        result_viper_wilcoxon$P.Value.I=two2one(result_viper_wilcoxon$P.Value, two = as.logical(result_viper_wilcoxon$two), invert = as.logical(result_viper_wilcoxon$invert_I))
        result_viper_wilcoxon_A=result_viper_wilcoxon[,c("geneset_id","P.Value.A")]
        result_viper_wilcoxon_I=result_viper_wilcoxon[,c("geneset_id","P.Value.I")]
        colnames(result_viper_wilcoxon_A)[2]=paste0(vector[i],".viper_wilcoxon")
        colnames(result_viper_wilcoxon_I)[2]=paste0(vector[i],".viper_wilcoxon")
        rm(result_viper_wilcoxon)
      }
      if("viper_kruskal" %in% enrichment_method){
        print("viper_kruskal will start.")
        result=Get_kruskal_test(geneset=geneset,group=group)
        result_viper_kruskal=result[,c("geneset_id","P.Value","logFC")]
        result_viper_kruskal$two=TRUE
        result_viper_kruskal$invert_A=ifelse(result_viper_kruskal$logFC<=0,TRUE,FALSE)
        result_viper_kruskal$invert_I=ifelse(result_viper_kruskal$logFC>0,TRUE,FALSE)
        result_viper_kruskal$P.Value.A=two2one(result_viper_kruskal$P.Value, two = as.logical(result_viper_kruskal$two), invert = as.logical(result_viper_kruskal$invert_A))
        result_viper_kruskal$P.Value.I=two2one(result_viper_kruskal$P.Value, two = as.logical(result_viper_kruskal$two), invert = as.logical(result_viper_kruskal$invert_I))
        result_viper_kruskal_A=result_viper_kruskal[,c("geneset_id","P.Value.A")]
        result_viper_kruskal_I=result_viper_kruskal[,c("geneset_id","P.Value.I")]
        colnames(result_viper_kruskal_A)[2]=paste0(vector[i],".viper_kruskal")
        colnames(result_viper_kruskal_I)[2]=paste0(vector[i],".viper_kruskal")
        rm(result_viper_kruskal)
      }
      if("viper_permutation" %in% enrichment_method){
        print("viper_permutation will start.")
        result=Get_permutation_test(geneset=geneset,group=group)
        result_viper_permutation=result[,c("geneset_id","P.Value","logFC")]
        result_viper_permutation$two=TRUE
        result_viper_permutation$invert_A=ifelse(result_viper_permutation$logFC<=0,TRUE,FALSE)
        result_viper_permutation$invert_I=ifelse(result_viper_permutation$logFC>0,TRUE,FALSE)
        result_viper_permutation=na.omit(result_viper_permutation)
        result_viper_permutation$P.Value.A=two2one(result_viper_permutation$P.Value, two = as.logical(result_viper_permutation$two), invert = as.logical(result_viper_permutation$invert_A))
        result_viper_permutation$P.Value.I=two2one(result_viper_permutation$P.Value, two = as.logical(result_viper_permutation$two), invert = as.logical(result_viper_permutation$invert_I))
        result_viper_permutation_A=result_viper_permutation[,c("geneset_id","P.Value.A")]
        result_viper_permutation_I=result_viper_permutation[,c("geneset_id","P.Value.I")]
        colnames(result_viper_permutation_A)[2]=paste0(vector[i],".viper_permutation")
        colnames(result_viper_permutation_I)[2]=paste0(vector[i],".viper_permutation")
        rm(result_viper_permutation)
      }
      if("viper_limma" %in% enrichment_method){
        print("viper_limma will start.")
        result=Get_limma(geneset=geneset,group=group)
        result_viper_limma=result[,c("geneset_id","P.Value","logFC")]
        result_viper_limma$two=TRUE
        result_viper_limma$invert_A=ifelse(result_viper_limma$logFC<=0,TRUE,FALSE)
        result_viper_limma$invert_I=ifelse(result_viper_limma$logFC>0,TRUE,FALSE)
        result_viper_limma$P.Value.A=two2one(result_viper_limma$P.Value, two = as.logical(result_viper_limma$two), invert = as.logical(result_viper_limma$invert_A))
        result_viper_limma$P.Value.I=two2one(result_viper_limma$P.Value, two = as.logical(result_viper_limma$two), invert = as.logical(result_viper_limma$invert_I))
        result_viper_limma_A=result_viper_limma[,c("geneset_id","P.Value.A")]
        result_viper_limma_I=result_viper_limma[,c("geneset_id","P.Value.I")]
        colnames(result_viper_limma_A)[2]=paste0(vector[i],".viper_limma")
        colnames(result_viper_limma_I)[2]=paste0(vector[i],".viper_limma")
        rm(result_viper_limma)
      }
      if("viper_anova" %in% enrichment_method){
        print("viper_anova will start.")
        result=Get_anova(geneset=geneset,group=group)
        result_viper_anova=result[,c("geneset_id","P.Value","logFC")]
        result_viper_anova$two=TRUE
        result_viper_anova$invert_A=ifelse(result_viper_anova$logFC<=0,TRUE,FALSE)
        result_viper_anova$invert_I=ifelse(result_viper_anova$logFC>0,TRUE,FALSE)
        result_viper_anova$P.Value.A=two2one(result_viper_anova$P.Value, two = as.logical(result_viper_anova$two), invert = as.logical(result_viper_anova$invert_A))
        result_viper_anova$P.Value.I=two2one(result_viper_anova$P.Value, two = as.logical(result_viper_anova$two), invert = as.logical(result_viper_anova$invert_I))
        result_viper_anova_A=result_viper_anova[,c("geneset_id","P.Value.A")]
        result_viper_anova_I=result_viper_anova[,c("geneset_id","P.Value.I")]
        colnames(result_viper_anova_A)[2]=paste0(vector[i],".viper_anova")
        colnames(result_viper_anova_I)[2]=paste0(vector[i],".viper_anova")
        rm(result_viper_anova)
      }
    } else {
      print("Not using viper.")
    }
    ##########################################wmean
    methods_to_check=c("wmean_t", "wmean_limma", "wmean_anova", "wmean_wilcoxon", "wmean_permutation", "wmean_kruskal")
    if (any(methods_to_check %in% enrichment_method)){
      geneSets_list=lapply(geneSets, function(x){ x@geneIds })
      names(geneSets_list)=names(geneSets)
      network=data.frame()
      for(exo in 1:length(geneSets_list)){
        geneSets_list_exo=as.data.frame(geneSets_list[[exo]])
        colnames(geneSets_list_exo)="target"
        geneSets_list_exo$source=names(geneSets_list)[exo]
        geneSets_list_exo=geneSets_list_exo[,c(2,1)]
        geneSets_list_exo$mor=1
        network=rbind(network,geneSets_list_exo)
      }
      network=network[network$target!="",]
      geneset=run_wmean(mat=exp,network=network,minsize = min.sz)
      geneset=geneset[geneset$statistic=="wmean",]
      geneset=geneset[,c("condition","source","score")]
      geneset=reshape2::dcast(geneset, source ~ condition, value.var = "score")
      rownames(geneset)=geneset[,1]
      geneset=geneset[,-1]
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("wmean_t" %in% enrichment_method){
        print("wmean_t will start.")
        result=Get_T_test(geneset=geneset,group=group)
        result_wmean_t=result[,c("geneset_id","P.Value","logFC")]
        result_wmean_t$two=TRUE
        result_wmean_t$invert_A=ifelse(result_wmean_t$logFC<=0,TRUE,FALSE)
        result_wmean_t$invert_I=ifelse(result_wmean_t$logFC>0,TRUE,FALSE)
        result_wmean_t$P.Value.A=two2one(result_wmean_t$P.Value, two = as.logical(result_wmean_t$two), invert = as.logical(result_wmean_t$invert_A))
        result_wmean_t$P.Value.I=two2one(result_wmean_t$P.Value, two = as.logical(result_wmean_t$two), invert = as.logical(result_wmean_t$invert_I))
        result_wmean_t_A=result_wmean_t[,c("geneset_id","P.Value.A")]
        result_wmean_t_I=result_wmean_t[,c("geneset_id","P.Value.I")]
        colnames(result_wmean_t_A)[2]=paste0(vector[i],".wmean_t")
        colnames(result_wmean_t_I)[2]=paste0(vector[i],".wmean_t")
        rm(result_wmean_t)
      }
      if("wmean_wilcoxon" %in% enrichment_method){
        print("wmean_wilcoxon will start.")
        result=Get_wilcox_test(geneset=geneset,group=group)
        result_wmean_wilcoxon=result[,c("geneset_id","P.Value","logFC")]
        result_wmean_wilcoxon$two=TRUE
        result_wmean_wilcoxon$invert_A=ifelse(result_wmean_wilcoxon$logFC<=0,TRUE,FALSE)
        result_wmean_wilcoxon$invert_I=ifelse(result_wmean_wilcoxon$logFC>0,TRUE,FALSE)
        result_wmean_wilcoxon$P.Value.A=two2one(result_wmean_wilcoxon$P.Value, two = as.logical(result_wmean_wilcoxon$two), invert = as.logical(result_wmean_wilcoxon$invert_A))
        result_wmean_wilcoxon$P.Value.I=two2one(result_wmean_wilcoxon$P.Value, two = as.logical(result_wmean_wilcoxon$two), invert = as.logical(result_wmean_wilcoxon$invert_I))
        result_wmean_wilcoxon_A=result_wmean_wilcoxon[,c("geneset_id","P.Value.A")]
        result_wmean_wilcoxon_I=result_wmean_wilcoxon[,c("geneset_id","P.Value.I")]
        colnames(result_wmean_wilcoxon_A)[2]=paste0(vector[i],".wmean_wilcoxon")
        colnames(result_wmean_wilcoxon_I)[2]=paste0(vector[i],".wmean_wilcoxon")
        rm(result_wmean_wilcoxon)
      }
      if("wmean_kruskal" %in% enrichment_method){
        print("wmean_kruskal will start.")
        result=Get_kruskal_test(geneset=geneset,group=group)
        result_wmean_kruskal=result[,c("geneset_id","P.Value","logFC")]
        result_wmean_kruskal$two=TRUE
        result_wmean_kruskal$invert_A=ifelse(result_wmean_kruskal$logFC<=0,TRUE,FALSE)
        result_wmean_kruskal$invert_I=ifelse(result_wmean_kruskal$logFC>0,TRUE,FALSE)
        result_wmean_kruskal$P.Value.A=two2one(result_wmean_kruskal$P.Value, two = as.logical(result_wmean_kruskal$two), invert = as.logical(result_wmean_kruskal$invert_A))
        result_wmean_kruskal$P.Value.I=two2one(result_wmean_kruskal$P.Value, two = as.logical(result_wmean_kruskal$two), invert = as.logical(result_wmean_kruskal$invert_I))
        result_wmean_kruskal_A=result_wmean_kruskal[,c("geneset_id","P.Value.A")]
        result_wmean_kruskal_I=result_wmean_kruskal[,c("geneset_id","P.Value.I")]
        colnames(result_wmean_kruskal_A)[2]=paste0(vector[i],".wmean_kruskal")
        colnames(result_wmean_kruskal_I)[2]=paste0(vector[i],".wmean_kruskal")
        rm(result_wmean_kruskal)
      }
      if("wmean_permutation" %in% enrichment_method){
        print("wmean_permutation will start.")
        result=Get_permutation_test(geneset=geneset,group=group)
        result_wmean_permutation=result[,c("geneset_id","P.Value","logFC")]
        result_wmean_permutation$two=TRUE
        result_wmean_permutation$invert_A=ifelse(result_wmean_permutation$logFC<=0,TRUE,FALSE)
        result_wmean_permutation$invert_I=ifelse(result_wmean_permutation$logFC>0,TRUE,FALSE)
        result_wmean_permutation=na.omit(result_wmean_permutation)
        result_wmean_permutation$P.Value.A=two2one(result_wmean_permutation$P.Value, two = as.logical(result_wmean_permutation$two), invert = as.logical(result_wmean_permutation$invert_A))
        result_wmean_permutation$P.Value.I=two2one(result_wmean_permutation$P.Value, two = as.logical(result_wmean_permutation$two), invert = as.logical(result_wmean_permutation$invert_I))
        result_wmean_permutation_A=result_wmean_permutation[,c("geneset_id","P.Value.A")]
        result_wmean_permutation_I=result_wmean_permutation[,c("geneset_id","P.Value.I")]
        colnames(result_wmean_permutation_A)[2]=paste0(vector[i],".wmean_permutation")
        colnames(result_wmean_permutation_I)[2]=paste0(vector[i],".wmean_permutation")
        rm(result_wmean_permutation)
      }
      if("wmean_limma" %in% enrichment_method){
        print("wmean_limma will start.")
        result=Get_limma(geneset=geneset,group=group)
        result_wmean_limma=result[,c("geneset_id","P.Value","logFC")]
        result_wmean_limma$two=TRUE
        result_wmean_limma$invert_A=ifelse(result_wmean_limma$logFC<=0,TRUE,FALSE)
        result_wmean_limma$invert_I=ifelse(result_wmean_limma$logFC>0,TRUE,FALSE)
        result_wmean_limma$P.Value.A=two2one(result_wmean_limma$P.Value, two = as.logical(result_wmean_limma$two), invert = as.logical(result_wmean_limma$invert_A))
        result_wmean_limma$P.Value.I=two2one(result_wmean_limma$P.Value, two = as.logical(result_wmean_limma$two), invert = as.logical(result_wmean_limma$invert_I))
        result_wmean_limma_A=result_wmean_limma[,c("geneset_id","P.Value.A")]
        result_wmean_limma_I=result_wmean_limma[,c("geneset_id","P.Value.I")]
        colnames(result_wmean_limma_A)[2]=paste0(vector[i],".wmean_limma")
        colnames(result_wmean_limma_I)[2]=paste0(vector[i],".wmean_limma")
        rm(result_wmean_limma)
      }
      if("wmean_anova" %in% enrichment_method){
        print("wmean_anova will start.")
        result=Get_anova(geneset=geneset,group=group)
        result_wmean_anova=result[,c("geneset_id","P.Value","logFC")]
        result_wmean_anova$two=TRUE
        result_wmean_anova$invert_A=ifelse(result_wmean_anova$logFC<=0,TRUE,FALSE)
        result_wmean_anova$invert_I=ifelse(result_wmean_anova$logFC>0,TRUE,FALSE)
        result_wmean_anova$P.Value.A=two2one(result_wmean_anova$P.Value, two = as.logical(result_wmean_anova$two), invert = as.logical(result_wmean_anova$invert_A))
        result_wmean_anova$P.Value.I=two2one(result_wmean_anova$P.Value, two = as.logical(result_wmean_anova$two), invert = as.logical(result_wmean_anova$invert_I))
        result_wmean_anova_A=result_wmean_anova[,c("geneset_id","P.Value.A")]
        result_wmean_anova_I=result_wmean_anova[,c("geneset_id","P.Value.I")]
        colnames(result_wmean_anova_A)[2]=paste0(vector[i],".wmean_anova")
        colnames(result_wmean_anova_I)[2]=paste0(vector[i],".wmean_anova")
        rm(result_wmean_anova)
      }
    } else {
      print("Not using wmean.")
    }
    ##########################################norm_wmean
    methods_to_check=c("norm_wmean_t", "norm_wmean_limma", "norm_wmean_anova", "norm_wmean_wilcoxon", "norm_wmean_permutation", "norm_wmean_kruskal")
    if (any(methods_to_check %in% enrichment_method)){
      geneSets_list=lapply(geneSets, function(x){ x@geneIds })
      names(geneSets_list)=names(geneSets)
      network=data.frame()
      for(exo in 1:length(geneSets_list)){
        geneSets_list_exo=as.data.frame(geneSets_list[[exo]])
        colnames(geneSets_list_exo)="target"
        geneSets_list_exo$source=names(geneSets_list)[exo]
        geneSets_list_exo=geneSets_list_exo[,c(2,1)]
        geneSets_list_exo$mor=1
        network=rbind(network,geneSets_list_exo)
      }
      network=network[network$target!="",]
      geneset=run_wmean(mat=exp,network=network,minsize = min.sz)
      geneset=geneset[geneset$statistic=="norm_wmean",]
      geneset=geneset[,c("condition","source","score")]
      geneset=reshape2::dcast(geneset, source ~ condition, value.var = "score")
      rownames(geneset)=geneset[,1]
      geneset=geneset[,-1]
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("norm_wmean_t" %in% enrichment_method){
        print("norm_wmean_t will start.")
        result=Get_T_test(geneset=geneset,group=group)
        result_norm_wmean_t=result[,c("geneset_id","P.Value","logFC")]
        result_norm_wmean_t$two=TRUE
        result_norm_wmean_t$invert_A=ifelse(result_norm_wmean_t$logFC<=0,TRUE,FALSE)
        result_norm_wmean_t$invert_I=ifelse(result_norm_wmean_t$logFC>0,TRUE,FALSE)
        result_norm_wmean_t$P.Value.A=two2one(result_norm_wmean_t$P.Value, two = as.logical(result_norm_wmean_t$two), invert = as.logical(result_norm_wmean_t$invert_A))
        result_norm_wmean_t$P.Value.I=two2one(result_norm_wmean_t$P.Value, two = as.logical(result_norm_wmean_t$two), invert = as.logical(result_norm_wmean_t$invert_I))
        result_norm_wmean_t_A=result_norm_wmean_t[,c("geneset_id","P.Value.A")]
        result_norm_wmean_t_I=result_norm_wmean_t[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wmean_t_A)[2]=paste0(vector[i],".norm_wmean_t")
        colnames(result_norm_wmean_t_I)[2]=paste0(vector[i],".norm_wmean_t")
        rm(result_norm_wmean_t)
      }
      if("norm_wmean_wilcoxon" %in% enrichment_method){
        print("norm_wmean_wilcoxon will start.")
        result=Get_wilcox_test(geneset=geneset,group=group)
        result_norm_wmean_wilcoxon=result[,c("geneset_id","P.Value","logFC")]
        result_norm_wmean_wilcoxon$two=TRUE
        result_norm_wmean_wilcoxon$invert_A=ifelse(result_norm_wmean_wilcoxon$logFC<=0,TRUE,FALSE)
        result_norm_wmean_wilcoxon$invert_I=ifelse(result_norm_wmean_wilcoxon$logFC>0,TRUE,FALSE)
        result_norm_wmean_wilcoxon$P.Value.A=two2one(result_norm_wmean_wilcoxon$P.Value, two = as.logical(result_norm_wmean_wilcoxon$two), invert = as.logical(result_norm_wmean_wilcoxon$invert_A))
        result_norm_wmean_wilcoxon$P.Value.I=two2one(result_norm_wmean_wilcoxon$P.Value, two = as.logical(result_norm_wmean_wilcoxon$two), invert = as.logical(result_norm_wmean_wilcoxon$invert_I))
        result_norm_wmean_wilcoxon_A=result_norm_wmean_wilcoxon[,c("geneset_id","P.Value.A")]
        result_norm_wmean_wilcoxon_I=result_norm_wmean_wilcoxon[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wmean_wilcoxon_A)[2]=paste0(vector[i],".norm_wmean_wilcoxon")
        colnames(result_norm_wmean_wilcoxon_I)[2]=paste0(vector[i],".norm_wmean_wilcoxon")
        rm(result_norm_wmean_wilcoxon)
      }
      if("norm_wmean_kruskal" %in% enrichment_method){
        print("norm_wmean_kruskal will start.")
        result=Get_kruskal_test(geneset=geneset,group=group)
        result_norm_wmean_kruskal=result[,c("geneset_id","P.Value","logFC")]
        result_norm_wmean_kruskal$two=TRUE
        result_norm_wmean_kruskal$invert_A=ifelse(result_norm_wmean_kruskal$logFC<=0,TRUE,FALSE)
        result_norm_wmean_kruskal$invert_I=ifelse(result_norm_wmean_kruskal$logFC>0,TRUE,FALSE)
        result_norm_wmean_kruskal$P.Value.A=two2one(result_norm_wmean_kruskal$P.Value, two = as.logical(result_norm_wmean_kruskal$two), invert = as.logical(result_norm_wmean_kruskal$invert_A))
        result_norm_wmean_kruskal$P.Value.I=two2one(result_norm_wmean_kruskal$P.Value, two = as.logical(result_norm_wmean_kruskal$two), invert = as.logical(result_norm_wmean_kruskal$invert_I))
        result_norm_wmean_kruskal_A=result_norm_wmean_kruskal[,c("geneset_id","P.Value.A")]
        result_norm_wmean_kruskal_I=result_norm_wmean_kruskal[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wmean_kruskal_A)[2]=paste0(vector[i],".norm_wmean_kruskal")
        colnames(result_norm_wmean_kruskal_I)[2]=paste0(vector[i],".norm_wmean_kruskal")
        rm(result_norm_wmean_kruskal)
      }
      if("norm_wmean_permutation" %in% enrichment_method){
        print("norm_wmean_permutation will start.")
        result=Get_permutation_test(geneset=geneset,group=group)
        result_norm_wmean_permutation=result[,c("geneset_id","P.Value","logFC")]
        result_norm_wmean_permutation$two=TRUE
        result_norm_wmean_permutation$invert_A=ifelse(result_norm_wmean_permutation$logFC<=0,TRUE,FALSE)
        result_norm_wmean_permutation$invert_I=ifelse(result_norm_wmean_permutation$logFC>0,TRUE,FALSE)
        result_norm_wmean_permutation=na.omit(result_norm_wmean_permutation)
        result_norm_wmean_permutation$P.Value.A=two2one(result_norm_wmean_permutation$P.Value, two = as.logical(result_norm_wmean_permutation$two), invert = as.logical(result_norm_wmean_permutation$invert_A))
        result_norm_wmean_permutation$P.Value.I=two2one(result_norm_wmean_permutation$P.Value, two = as.logical(result_norm_wmean_permutation$two), invert = as.logical(result_norm_wmean_permutation$invert_I))
        result_norm_wmean_permutation_A=result_norm_wmean_permutation[,c("geneset_id","P.Value.A")]
        result_norm_wmean_permutation_I=result_norm_wmean_permutation[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wmean_permutation_A)[2]=paste0(vector[i],".norm_wmean_permutation")
        colnames(result_norm_wmean_permutation_I)[2]=paste0(vector[i],".norm_wmean_permutation")
        rm(result_norm_wmean_permutation)
      }
      if("norm_wmean_limma" %in% enrichment_method){
        print("norm_wmean_limma will start.")
        result=Get_limma(geneset=geneset,group=group)
        result_norm_wmean_limma=result[,c("geneset_id","P.Value","logFC")]
        result_norm_wmean_limma$two=TRUE
        result_norm_wmean_limma$invert_A=ifelse(result_norm_wmean_limma$logFC<=0,TRUE,FALSE)
        result_norm_wmean_limma$invert_I=ifelse(result_norm_wmean_limma$logFC>0,TRUE,FALSE)
        result_norm_wmean_limma$P.Value.A=two2one(result_norm_wmean_limma$P.Value, two = as.logical(result_norm_wmean_limma$two), invert = as.logical(result_norm_wmean_limma$invert_A))
        result_norm_wmean_limma$P.Value.I=two2one(result_norm_wmean_limma$P.Value, two = as.logical(result_norm_wmean_limma$two), invert = as.logical(result_norm_wmean_limma$invert_I))
        result_norm_wmean_limma_A=result_norm_wmean_limma[,c("geneset_id","P.Value.A")]
        result_norm_wmean_limma_I=result_norm_wmean_limma[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wmean_limma_A)[2]=paste0(vector[i],".norm_wmean_limma")
        colnames(result_norm_wmean_limma_I)[2]=paste0(vector[i],".norm_wmean_limma")
        rm(result_norm_wmean_limma)
      }
      if("norm_wmean_anova" %in% enrichment_method){
        print("norm_wmean_anova will start.")
        result=Get_anova(geneset=geneset,group=group)
        result_norm_wmean_anova=result[,c("geneset_id","P.Value","logFC")]
        result_norm_wmean_anova$two=TRUE
        result_norm_wmean_anova$invert_A=ifelse(result_norm_wmean_anova$logFC<=0,TRUE,FALSE)
        result_norm_wmean_anova$invert_I=ifelse(result_norm_wmean_anova$logFC>0,TRUE,FALSE)
        result_norm_wmean_anova$P.Value.A=two2one(result_norm_wmean_anova$P.Value, two = as.logical(result_norm_wmean_anova$two), invert = as.logical(result_norm_wmean_anova$invert_A))
        result_norm_wmean_anova$P.Value.I=two2one(result_norm_wmean_anova$P.Value, two = as.logical(result_norm_wmean_anova$two), invert = as.logical(result_norm_wmean_anova$invert_I))
        result_norm_wmean_anova_A=result_norm_wmean_anova[,c("geneset_id","P.Value.A")]
        result_norm_wmean_anova_I=result_norm_wmean_anova[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wmean_anova_A)[2]=paste0(vector[i],".norm_wmean_anova")
        colnames(result_norm_wmean_anova_I)[2]=paste0(vector[i],".norm_wmean_anova")
        rm(result_norm_wmean_anova)
      }
    } else {
      print("Not using norm_wmean.")
    }
    ##########################################corr_wmean
    methods_to_check=c("corr_wmean_t", "corr_wmean_limma", "corr_wmean_anova", "corr_wmean_wilcoxon", "corr_wmean_permutation", "corr_wmean_kruskal")
    if (any(methods_to_check %in% enrichment_method)){
      geneSets_list=lapply(geneSets, function(x){ x@geneIds })
      names(geneSets_list)=names(geneSets)
      network=data.frame()
      for(exo in 1:length(geneSets_list)){
        geneSets_list_exo=as.data.frame(geneSets_list[[exo]])
        colnames(geneSets_list_exo)="target"
        geneSets_list_exo$source=names(geneSets_list)[exo]
        geneSets_list_exo=geneSets_list_exo[,c(2,1)]
        geneSets_list_exo$mor=1
        network=rbind(network,geneSets_list_exo)
      }
      network=network[network$target!="",]
      geneset=run_wmean(mat=exp,network=network,minsize = min.sz)
      geneset=geneset[geneset$statistic=="corr_wmean",]
      geneset=geneset[,c("condition","source","score")]
      geneset=reshape2::dcast(geneset, source ~ condition, value.var = "score")
      rownames(geneset)=geneset[,1]
      geneset=geneset[,-1]
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("corr_wmean_t" %in% enrichment_method){
        print("corr_wmean_t will start.")
        result=Get_T_test(geneset=geneset,group=group)
        result_corr_wmean_t=result[,c("geneset_id","P.Value","logFC")]
        result_corr_wmean_t$two=TRUE
        result_corr_wmean_t$invert_A=ifelse(result_corr_wmean_t$logFC<=0,TRUE,FALSE)
        result_corr_wmean_t$invert_I=ifelse(result_corr_wmean_t$logFC>0,TRUE,FALSE)
        result_corr_wmean_t$P.Value.A=two2one(result_corr_wmean_t$P.Value, two = as.logical(result_corr_wmean_t$two), invert = as.logical(result_corr_wmean_t$invert_A))
        result_corr_wmean_t$P.Value.I=two2one(result_corr_wmean_t$P.Value, two = as.logical(result_corr_wmean_t$two), invert = as.logical(result_corr_wmean_t$invert_I))
        result_corr_wmean_t_A=result_corr_wmean_t[,c("geneset_id","P.Value.A")]
        result_corr_wmean_t_I=result_corr_wmean_t[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wmean_t_A)[2]=paste0(vector[i],".corr_wmean_t")
        colnames(result_corr_wmean_t_I)[2]=paste0(vector[i],".corr_wmean_t")
        rm(result_corr_wmean_t)
      }
      if("corr_wmean_wilcoxon" %in% enrichment_method){
        print("corr_wmean_wilcoxon will start.")
        result=Get_wilcox_test(geneset=geneset,group=group)
        result_corr_wmean_wilcoxon=result[,c("geneset_id","P.Value","logFC")]
        result_corr_wmean_wilcoxon$two=TRUE
        result_corr_wmean_wilcoxon$invert_A=ifelse(result_corr_wmean_wilcoxon$logFC<=0,TRUE,FALSE)
        result_corr_wmean_wilcoxon$invert_I=ifelse(result_corr_wmean_wilcoxon$logFC>0,TRUE,FALSE)
        result_corr_wmean_wilcoxon$P.Value.A=two2one(result_corr_wmean_wilcoxon$P.Value, two = as.logical(result_corr_wmean_wilcoxon$two), invert = as.logical(result_corr_wmean_wilcoxon$invert_A))
        result_corr_wmean_wilcoxon$P.Value.I=two2one(result_corr_wmean_wilcoxon$P.Value, two = as.logical(result_corr_wmean_wilcoxon$two), invert = as.logical(result_corr_wmean_wilcoxon$invert_I))
        result_corr_wmean_wilcoxon_A=result_corr_wmean_wilcoxon[,c("geneset_id","P.Value.A")]
        result_corr_wmean_wilcoxon_I=result_corr_wmean_wilcoxon[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wmean_wilcoxon_A)[2]=paste0(vector[i],".corr_wmean_wilcoxon")
        colnames(result_corr_wmean_wilcoxon_I)[2]=paste0(vector[i],".corr_wmean_wilcoxon")
        rm(result_corr_wmean_wilcoxon)
      }
      if("corr_wmean_kruskal" %in% enrichment_method){
        print("corr_wmean_kruskal will start.")
        result=Get_kruskal_test(geneset=geneset,group=group)
        result_corr_wmean_kruskal=result[,c("geneset_id","P.Value","logFC")]
        result_corr_wmean_kruskal$two=TRUE
        result_corr_wmean_kruskal$invert_A=ifelse(result_corr_wmean_kruskal$logFC<=0,TRUE,FALSE)
        result_corr_wmean_kruskal$invert_I=ifelse(result_corr_wmean_kruskal$logFC>0,TRUE,FALSE)
        result_corr_wmean_kruskal$P.Value.A=two2one(result_corr_wmean_kruskal$P.Value, two = as.logical(result_corr_wmean_kruskal$two), invert = as.logical(result_corr_wmean_kruskal$invert_A))
        result_corr_wmean_kruskal$P.Value.I=two2one(result_corr_wmean_kruskal$P.Value, two = as.logical(result_corr_wmean_kruskal$two), invert = as.logical(result_corr_wmean_kruskal$invert_I))
        result_corr_wmean_kruskal_A=result_corr_wmean_kruskal[,c("geneset_id","P.Value.A")]
        result_corr_wmean_kruskal_I=result_corr_wmean_kruskal[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wmean_kruskal_A)[2]=paste0(vector[i],".corr_wmean_kruskal")
        colnames(result_corr_wmean_kruskal_I)[2]=paste0(vector[i],".corr_wmean_kruskal")
        rm(result_corr_wmean_kruskal)
      }
      if("corr_wmean_permutation" %in% enrichment_method){
        print("corr_wmean_permutation will start.")
        result=Get_permutation_test(geneset=geneset,group=group)
        result_corr_wmean_permutation=result[,c("geneset_id","P.Value","logFC")]
        result_corr_wmean_permutation$two=TRUE
        result_corr_wmean_permutation$invert_A=ifelse(result_corr_wmean_permutation$logFC<=0,TRUE,FALSE)
        result_corr_wmean_permutation$invert_I=ifelse(result_corr_wmean_permutation$logFC>0,TRUE,FALSE)
        result_corr_wmean_permutation=na.omit(result_corr_wmean_permutation)
        result_corr_wmean_permutation$P.Value.A=two2one(result_corr_wmean_permutation$P.Value, two = as.logical(result_corr_wmean_permutation$two), invert = as.logical(result_corr_wmean_permutation$invert_A))
        result_corr_wmean_permutation$P.Value.I=two2one(result_corr_wmean_permutation$P.Value, two = as.logical(result_corr_wmean_permutation$two), invert = as.logical(result_corr_wmean_permutation$invert_I))
        result_corr_wmean_permutation_A=result_corr_wmean_permutation[,c("geneset_id","P.Value.A")]
        result_corr_wmean_permutation_I=result_corr_wmean_permutation[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wmean_permutation_A)[2]=paste0(vector[i],".corr_wmean_permutation")
        colnames(result_corr_wmean_permutation_I)[2]=paste0(vector[i],".corr_wmean_permutation")
        rm(result_corr_wmean_permutation)
      }
      if("corr_wmean_limma" %in% enrichment_method){
        print("corr_wmean_limma will start.")
        result=Get_limma(geneset=geneset,group=group)
        result_corr_wmean_limma=result[,c("geneset_id","P.Value","logFC")]
        result_corr_wmean_limma$two=TRUE
        result_corr_wmean_limma$invert_A=ifelse(result_corr_wmean_limma$logFC<=0,TRUE,FALSE)
        result_corr_wmean_limma$invert_I=ifelse(result_corr_wmean_limma$logFC>0,TRUE,FALSE)
        result_corr_wmean_limma$P.Value.A=two2one(result_corr_wmean_limma$P.Value, two = as.logical(result_corr_wmean_limma$two), invert = as.logical(result_corr_wmean_limma$invert_A))
        result_corr_wmean_limma$P.Value.I=two2one(result_corr_wmean_limma$P.Value, two = as.logical(result_corr_wmean_limma$two), invert = as.logical(result_corr_wmean_limma$invert_I))
        result_corr_wmean_limma_A=result_corr_wmean_limma[,c("geneset_id","P.Value.A")]
        result_corr_wmean_limma_I=result_corr_wmean_limma[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wmean_limma_A)[2]=paste0(vector[i],".corr_wmean_limma")
        colnames(result_corr_wmean_limma_I)[2]=paste0(vector[i],".corr_wmean_limma")
        rm(result_corr_wmean_limma)
      }
      if("corr_wmean_anova" %in% enrichment_method){
        print("corr_wmean_anova will start.")
        result=Get_anova(geneset=geneset,group=group)
        result_corr_wmean_anova=result[,c("geneset_id","P.Value","logFC")]
        result_corr_wmean_anova$two=TRUE
        result_corr_wmean_anova$invert_A=ifelse(result_corr_wmean_anova$logFC<=0,TRUE,FALSE)
        result_corr_wmean_anova$invert_I=ifelse(result_corr_wmean_anova$logFC>0,TRUE,FALSE)
        result_corr_wmean_anova$P.Value.A=two2one(result_corr_wmean_anova$P.Value, two = as.logical(result_corr_wmean_anova$two), invert = as.logical(result_corr_wmean_anova$invert_A))
        result_corr_wmean_anova$P.Value.I=two2one(result_corr_wmean_anova$P.Value, two = as.logical(result_corr_wmean_anova$two), invert = as.logical(result_corr_wmean_anova$invert_I))
        result_corr_wmean_anova_A=result_corr_wmean_anova[,c("geneset_id","P.Value.A")]
        result_corr_wmean_anova_I=result_corr_wmean_anova[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wmean_anova_A)[2]=paste0(vector[i],".corr_wmean_anova")
        colnames(result_corr_wmean_anova_I)[2]=paste0(vector[i],".corr_wmean_anova")
        rm(result_corr_wmean_anova)
      }
    } else {
      print("Not using corr_wmean.")
    }
    ##########################################wsum
    methods_to_check=c("wsum_t", "wsum_limma", "wsum_anova", "wsum_wilcoxon", "wsum_permutation", "wsum_kruskal")
    if (any(methods_to_check %in% enrichment_method)){
      geneSets_list=lapply(geneSets, function(x){ x@geneIds })
      names(geneSets_list)=names(geneSets)
      network=data.frame()
      for(exo in 1:length(geneSets_list)){
        geneSets_list_exo=as.data.frame(geneSets_list[[exo]])
        colnames(geneSets_list_exo)="target"
        geneSets_list_exo$source=names(geneSets_list)[exo]
        geneSets_list_exo=geneSets_list_exo[,c(2,1)]
        geneSets_list_exo$mor=1
        network=rbind(network,geneSets_list_exo)
      }
      network=network[network$target!="",]
      geneset=run_wsum(mat=exp,network=network,minsize = min.sz)
      geneset=geneset[geneset$statistic=="wsum",]
      geneset=geneset[,c("condition","source","score")]
      geneset=reshape2::dcast(geneset, source ~ condition, value.var = "score")
      rownames(geneset)=geneset[,1]
      geneset=geneset[,-1]
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("wsum_t" %in% enrichment_method){
        print("wsum_t will start.")
        result=Get_T_test(geneset=geneset,group=group)
        result_wsum_t=result[,c("geneset_id","P.Value","logFC")]
        result_wsum_t$two=TRUE
        result_wsum_t$invert_A=ifelse(result_wsum_t$logFC<=0,TRUE,FALSE)
        result_wsum_t$invert_I=ifelse(result_wsum_t$logFC>0,TRUE,FALSE)
        result_wsum_t$P.Value.A=two2one(result_wsum_t$P.Value, two = as.logical(result_wsum_t$two), invert = as.logical(result_wsum_t$invert_A))
        result_wsum_t$P.Value.I=two2one(result_wsum_t$P.Value, two = as.logical(result_wsum_t$two), invert = as.logical(result_wsum_t$invert_I))
        result_wsum_t_A=result_wsum_t[,c("geneset_id","P.Value.A")]
        result_wsum_t_I=result_wsum_t[,c("geneset_id","P.Value.I")]
        colnames(result_wsum_t_A)[2]=paste0(vector[i],".wsum_t")
        colnames(result_wsum_t_I)[2]=paste0(vector[i],".wsum_t")
        rm(result_wsum_t)
      }
      if("wsum_wilcoxon" %in% enrichment_method){
        print("wsum_wilcoxon will start.")
        result=Get_wilcox_test(geneset=geneset,group=group)
        result_wsum_wilcoxon=result[,c("geneset_id","P.Value","logFC")]
        result_wsum_wilcoxon$two=TRUE
        result_wsum_wilcoxon$invert_A=ifelse(result_wsum_wilcoxon$logFC<=0,TRUE,FALSE)
        result_wsum_wilcoxon$invert_I=ifelse(result_wsum_wilcoxon$logFC>0,TRUE,FALSE)
        result_wsum_wilcoxon$P.Value.A=two2one(result_wsum_wilcoxon$P.Value, two = as.logical(result_wsum_wilcoxon$two), invert = as.logical(result_wsum_wilcoxon$invert_A))
        result_wsum_wilcoxon$P.Value.I=two2one(result_wsum_wilcoxon$P.Value, two = as.logical(result_wsum_wilcoxon$two), invert = as.logical(result_wsum_wilcoxon$invert_I))
        result_wsum_wilcoxon_A=result_wsum_wilcoxon[,c("geneset_id","P.Value.A")]
        result_wsum_wilcoxon_I=result_wsum_wilcoxon[,c("geneset_id","P.Value.I")]
        colnames(result_wsum_wilcoxon_A)[2]=paste0(vector[i],".wsum_wilcoxon")
        colnames(result_wsum_wilcoxon_I)[2]=paste0(vector[i],".wsum_wilcoxon")
        rm(result_wsum_wilcoxon)
      }
      if("wsum_kruskal" %in% enrichment_method){
        print("wsum_kruskal will start.")
        result=Get_kruskal_test(geneset=geneset,group=group)
        result_wsum_kruskal=result[,c("geneset_id","P.Value","logFC")]
        result_wsum_kruskal$two=TRUE
        result_wsum_kruskal$invert_A=ifelse(result_wsum_kruskal$logFC<=0,TRUE,FALSE)
        result_wsum_kruskal$invert_I=ifelse(result_wsum_kruskal$logFC>0,TRUE,FALSE)
        result_wsum_kruskal$P.Value.A=two2one(result_wsum_kruskal$P.Value, two = as.logical(result_wsum_kruskal$two), invert = as.logical(result_wsum_kruskal$invert_A))
        result_wsum_kruskal$P.Value.I=two2one(result_wsum_kruskal$P.Value, two = as.logical(result_wsum_kruskal$two), invert = as.logical(result_wsum_kruskal$invert_I))
        result_wsum_kruskal_A=result_wsum_kruskal[,c("geneset_id","P.Value.A")]
        result_wsum_kruskal_I=result_wsum_kruskal[,c("geneset_id","P.Value.I")]
        colnames(result_wsum_kruskal_A)[2]=paste0(vector[i],".wsum_kruskal")
        colnames(result_wsum_kruskal_I)[2]=paste0(vector[i],".wsum_kruskal")
        rm(result_wsum_kruskal)
      }
      if("wsum_permutation" %in% enrichment_method){
        print("wsum_permutation will start.")
        result=Get_permutation_test(geneset=geneset,group=group)
        result_wsum_permutation=result[,c("geneset_id","P.Value","logFC")]
        result_wsum_permutation$two=TRUE
        result_wsum_permutation$invert_A=ifelse(result_wsum_permutation$logFC<=0,TRUE,FALSE)
        result_wsum_permutation$invert_I=ifelse(result_wsum_permutation$logFC>0,TRUE,FALSE)
        result_wsum_permutation=na.omit(result_wsum_permutation)
        result_wsum_permutation$P.Value.A=two2one(result_wsum_permutation$P.Value, two = as.logical(result_wsum_permutation$two), invert = as.logical(result_wsum_permutation$invert_A))
        result_wsum_permutation$P.Value.I=two2one(result_wsum_permutation$P.Value, two = as.logical(result_wsum_permutation$two), invert = as.logical(result_wsum_permutation$invert_I))
        result_wsum_permutation_A=result_wsum_permutation[,c("geneset_id","P.Value.A")]
        result_wsum_permutation_I=result_wsum_permutation[,c("geneset_id","P.Value.I")]
        colnames(result_wsum_permutation_A)[2]=paste0(vector[i],".wsum_permutation")
        colnames(result_wsum_permutation_I)[2]=paste0(vector[i],".wsum_permutation")
        rm(result_wsum_kruskal)
      }
      if("wsum_limma" %in% enrichment_method){
        print("wsum_limma will start.")
        result=Get_limma(geneset=geneset,group=group)
        result_wsum_limma=result[,c("geneset_id","P.Value","logFC")]
        result_wsum_limma$two=TRUE
        result_wsum_limma$invert_A=ifelse(result_wsum_limma$logFC<=0,TRUE,FALSE)
        result_wsum_limma$invert_I=ifelse(result_wsum_limma$logFC>0,TRUE,FALSE)
        result_wsum_limma$P.Value.A=two2one(result_wsum_limma$P.Value, two = as.logical(result_wsum_limma$two), invert = as.logical(result_wsum_limma$invert_A))
        result_wsum_limma$P.Value.I=two2one(result_wsum_limma$P.Value, two = as.logical(result_wsum_limma$two), invert = as.logical(result_wsum_limma$invert_I))
        result_wsum_limma_A=result_wsum_limma[,c("geneset_id","P.Value.A")]
        result_wsum_limma_I=result_wsum_limma[,c("geneset_id","P.Value.I")]
        colnames(result_wsum_limma_A)[2]=paste0(vector[i],".wsum_limma")
        colnames(result_wsum_limma_I)[2]=paste0(vector[i],".wsum_limma")
        rm(result_wsum_limma)
      }
      if("wsum_anova" %in% enrichment_method){
        print("wsum_anova will start.")
        result=Get_anova(geneset=geneset,group=group)
        result_wsum_anova=result[,c("geneset_id","P.Value","logFC")]
        result_wsum_anova$two=TRUE
        result_wsum_anova$invert_A=ifelse(result_wsum_anova$logFC<=0,TRUE,FALSE)
        result_wsum_anova$invert_I=ifelse(result_wsum_anova$logFC>0,TRUE,FALSE)
        result_wsum_anova$P.Value.A=two2one(result_wsum_anova$P.Value, two = as.logical(result_wsum_anova$two), invert = as.logical(result_wsum_anova$invert_A))
        result_wsum_anova$P.Value.I=two2one(result_wsum_anova$P.Value, two = as.logical(result_wsum_anova$two), invert = as.logical(result_wsum_anova$invert_I))
        result_wsum_anova_A=result_wsum_anova[,c("geneset_id","P.Value.A")]
        result_wsum_anova_I=result_wsum_anova[,c("geneset_id","P.Value.I")]
        colnames(result_wsum_anova_A)[2]=paste0(vector[i],".wsum_anova")
        colnames(result_wsum_anova_I)[2]=paste0(vector[i],".wsum_anova")
        rm(result_wsum_anova)
      }
    } else {
      print("Not using wsum.")
    }
    ##########################################norm_wsum
    methods_to_check=c("norm_wsum_t", "norm_wsum_limma", "norm_wsum_anova", "norm_wsum_wilcoxon", "norm_wsum_permutation", "norm_wsum_kruskal")
    if (any(methods_to_check %in% enrichment_method)){
      geneSets_list=lapply(geneSets, function(x){ x@geneIds })
      names(geneSets_list)=names(geneSets)
      network=data.frame()
      for(exo in 1:length(geneSets_list)){
        geneSets_list_exo=as.data.frame(geneSets_list[[exo]])
        colnames(geneSets_list_exo)="target"
        geneSets_list_exo$source=names(geneSets_list)[exo]
        geneSets_list_exo=geneSets_list_exo[,c(2,1)]
        geneSets_list_exo$mor=1
        network=rbind(network,geneSets_list_exo)
      }
      network=network[network$target!="",]
      geneset=run_wsum(mat=exp,network=network,minsize = min.sz)
      geneset=geneset[geneset$statistic=="norm_wsum",]
      geneset=geneset[,c("condition","source","score")]
      geneset=reshape2::dcast(geneset, source ~ condition, value.var = "score")
      rownames(geneset)=geneset[,1]
      geneset=geneset[,-1]
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("norm_wsum_t" %in% enrichment_method){
        print("norm_wsum_t will start.")
        result=Get_T_test(geneset=geneset,group=group)
        result_norm_wsum_t=result[,c("geneset_id","P.Value","logFC")]
        result_norm_wsum_t$two=TRUE
        result_norm_wsum_t$invert_A=ifelse(result_norm_wsum_t$logFC<=0,TRUE,FALSE)
        result_norm_wsum_t$invert_I=ifelse(result_norm_wsum_t$logFC>0,TRUE,FALSE)
        result_norm_wsum_t$P.Value.A=two2one(result_norm_wsum_t$P.Value, two = as.logical(result_norm_wsum_t$two), invert = as.logical(result_norm_wsum_t$invert_A))
        result_norm_wsum_t$P.Value.I=two2one(result_norm_wsum_t$P.Value, two = as.logical(result_norm_wsum_t$two), invert = as.logical(result_norm_wsum_t$invert_I))
        result_norm_wsum_t_A=result_norm_wsum_t[,c("geneset_id","P.Value.A")]
        result_norm_wsum_t_I=result_norm_wsum_t[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wsum_t_A)[2]=paste0(vector[i],".norm_wsum_t")
        colnames(result_norm_wsum_t_I)[2]=paste0(vector[i],".norm_wsum_t")
        rm(result_norm_wsum_t)
      }
      if("norm_wsum_wilcoxon" %in% enrichment_method){
        print("norm_wsum_wilcoxon will start.")
        result=Get_wilcox_test(geneset=geneset,group=group)
        result_norm_wsum_wilcoxon=result[,c("geneset_id","P.Value","logFC")]
        result_norm_wsum_wilcoxon$two=TRUE
        result_norm_wsum_wilcoxon$invert_A=ifelse(result_norm_wsum_wilcoxon$logFC<=0,TRUE,FALSE)
        result_norm_wsum_wilcoxon$invert_I=ifelse(result_norm_wsum_wilcoxon$logFC>0,TRUE,FALSE)
        result_norm_wsum_wilcoxon$P.Value.A=two2one(result_norm_wsum_wilcoxon$P.Value, two = as.logical(result_norm_wsum_wilcoxon$two), invert = as.logical(result_norm_wsum_wilcoxon$invert_A))
        result_norm_wsum_wilcoxon$P.Value.I=two2one(result_norm_wsum_wilcoxon$P.Value, two = as.logical(result_norm_wsum_wilcoxon$two), invert = as.logical(result_norm_wsum_wilcoxon$invert_I))
        result_norm_wsum_wilcoxon_A=result_norm_wsum_wilcoxon[,c("geneset_id","P.Value.A")]
        result_norm_wsum_wilcoxon_I=result_norm_wsum_wilcoxon[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wsum_wilcoxon_A)[2]=paste0(vector[i],".norm_wsum_wilcoxon")
        colnames(result_norm_wsum_wilcoxon_I)[2]=paste0(vector[i],".norm_wsum_wilcoxon")
        rm(result_norm_wsum_wilcoxon)
      }
      if("norm_wsum_kruskal" %in% enrichment_method){
        print("norm_wsum_kruskal will start.")
        result=Get_kruskal_test(geneset=geneset,group=group)
        result_norm_wsum_kruskal=result[,c("geneset_id","P.Value","logFC")]
        result_norm_wsum_kruskal$two=TRUE
        result_norm_wsum_kruskal$invert_A=ifelse(result_norm_wsum_kruskal$logFC<=0,TRUE,FALSE)
        result_norm_wsum_kruskal$invert_I=ifelse(result_norm_wsum_kruskal$logFC>0,TRUE,FALSE)
        result_norm_wsum_kruskal$P.Value.A=two2one(result_norm_wsum_kruskal$P.Value, two = as.logical(result_norm_wsum_kruskal$two), invert = as.logical(result_norm_wsum_kruskal$invert_A))
        result_norm_wsum_kruskal$P.Value.I=two2one(result_norm_wsum_kruskal$P.Value, two = as.logical(result_norm_wsum_kruskal$two), invert = as.logical(result_norm_wsum_kruskal$invert_I))
        result_norm_wsum_kruskal_A=result_norm_wsum_kruskal[,c("geneset_id","P.Value.A")]
        result_norm_wsum_kruskal_I=result_norm_wsum_kruskal[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wsum_kruskal_A)[2]=paste0(vector[i],".norm_wsum_kruskal")
        colnames(result_norm_wsum_kruskal_I)[2]=paste0(vector[i],".norm_wsum_kruskal")
        rm(result_norm_wsum_kruskal)
      }
      if("norm_wsum_permutation" %in% enrichment_method){
        print("norm_wsum_permutation will start.")
        result=Get_permutation_test(geneset=geneset,group=group)
        result_norm_wsum_permutation=result[,c("geneset_id","P.Value","logFC")]
        result_norm_wsum_permutation$two=TRUE
        result_norm_wsum_permutation$invert_A=ifelse(result_norm_wsum_permutation$logFC<=0,TRUE,FALSE)
        result_norm_wsum_permutation$invert_I=ifelse(result_norm_wsum_permutation$logFC>0,TRUE,FALSE)
        result_norm_wsum_permutation=na.omit(result_norm_wsum_permutation)
        result_norm_wsum_permutation$P.Value.A=two2one(result_norm_wsum_permutation$P.Value, two = as.logical(result_norm_wsum_permutation$two), invert = as.logical(result_norm_wsum_permutation$invert_A))
        result_norm_wsum_permutation$P.Value.I=two2one(result_norm_wsum_permutation$P.Value, two = as.logical(result_norm_wsum_permutation$two), invert = as.logical(result_norm_wsum_permutation$invert_I))
        result_norm_wsum_permutation_A=result_norm_wsum_permutation[,c("geneset_id","P.Value.A")]
        result_norm_wsum_permutation_I=result_norm_wsum_permutation[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wsum_permutation_A)[2]=paste0(vector[i],".norm_wsum_permutation")
        colnames(result_norm_wsum_permutation_I)[2]=paste0(vector[i],".norm_wsum_permutation")
        rm(result_norm_wsum_permutation)
      }
      if("norm_wsum_limma" %in% enrichment_method){
        print("norm_wsum_limma will start.")
        result=Get_limma(geneset=geneset,group=group)
        result_norm_wsum_limma=result[,c("geneset_id","P.Value","logFC")]
        result_norm_wsum_limma$two=TRUE
        result_norm_wsum_limma$invert_A=ifelse(result_norm_wsum_limma$logFC<=0,TRUE,FALSE)
        result_norm_wsum_limma$invert_I=ifelse(result_norm_wsum_limma$logFC>0,TRUE,FALSE)
        result_norm_wsum_limma$P.Value.A=two2one(result_norm_wsum_limma$P.Value, two = as.logical(result_norm_wsum_limma$two), invert = as.logical(result_norm_wsum_limma$invert_A))
        result_norm_wsum_limma$P.Value.I=two2one(result_norm_wsum_limma$P.Value, two = as.logical(result_norm_wsum_limma$two), invert = as.logical(result_norm_wsum_limma$invert_I))
        result_norm_wsum_limma_A=result_norm_wsum_limma[,c("geneset_id","P.Value.A")]
        result_norm_wsum_limma_I=result_norm_wsum_limma[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wsum_limma_A)[2]=paste0(vector[i],".norm_wsum_limma")
        colnames(result_norm_wsum_limma_I)[2]=paste0(vector[i],".norm_wsum_limma")
        rm(result_norm_wsum_limma)
      }
      if("norm_wsum_anova" %in% enrichment_method){
        print("norm_wsum_anova will start.")
        result=Get_anova(geneset=geneset,group=group)
        result_norm_wsum_anova=result[,c("geneset_id","P.Value","logFC")]
        result_norm_wsum_anova$two=TRUE
        result_norm_wsum_anova$invert_A=ifelse(result_norm_wsum_anova$logFC<=0,TRUE,FALSE)
        result_norm_wsum_anova$invert_I=ifelse(result_norm_wsum_anova$logFC>0,TRUE,FALSE)
        result_norm_wsum_anova$P.Value.A=two2one(result_norm_wsum_anova$P.Value, two = as.logical(result_norm_wsum_anova$two), invert = as.logical(result_norm_wsum_anova$invert_A))
        result_norm_wsum_anova$P.Value.I=two2one(result_norm_wsum_anova$P.Value, two = as.logical(result_norm_wsum_anova$two), invert = as.logical(result_norm_wsum_anova$invert_I))
        result_norm_wsum_anova_A=result_norm_wsum_anova[,c("geneset_id","P.Value.A")]
        result_norm_wsum_anova_I=result_norm_wsum_anova[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wsum_anova_A)[2]=paste0(vector[i],".norm_wsum_anova")
        colnames(result_norm_wsum_anova_I)[2]=paste0(vector[i],".norm_wsum_anova")
        rm(result_norm_wsum_anova)
      }
    } else {
      print("Not using norm_wsum.")
    }
    ##########################################corr_wsum
    methods_to_check=c("corr_wsum_t", "corr_wsum_limma", "corr_wsum_anova", "corr_wsum_wilcoxon", "corr_wsum_permutation", "corr_wsum_kruskal")
    if (any(methods_to_check %in% enrichment_method)){
      geneSets_list=lapply(geneSets, function(x){ x@geneIds })
      names(geneSets_list)=names(geneSets)
      network=data.frame()
      for(exo in 1:length(geneSets_list)){
        geneSets_list_exo=as.data.frame(geneSets_list[[exo]])
        colnames(geneSets_list_exo)="target"
        geneSets_list_exo$source=names(geneSets_list)[exo]
        geneSets_list_exo=geneSets_list_exo[,c(2,1)]
        geneSets_list_exo$mor=1
        network=rbind(network,geneSets_list_exo)
      }
      network=network[network$target!="",]
      geneset=run_wsum(mat=exp,network=network,minsize = min.sz)
      geneset=geneset[geneset$statistic=="corr_wsum",]
      geneset=geneset[,c("condition","source","score")]
      geneset=reshape2::dcast(geneset, source ~ condition, value.var = "score")
      rownames(geneset)=geneset[,1]
      geneset=geneset[,-1]
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("corr_wsum_t" %in% enrichment_method){
        print("corr_wsum_t will start.")
        result=Get_T_test(geneset=geneset,group=group)
        result_corr_wsum_t=result[,c("geneset_id","P.Value","logFC")]
        result_corr_wsum_t$two=TRUE
        result_corr_wsum_t$invert_A=ifelse(result_corr_wsum_t$logFC<=0,TRUE,FALSE)
        result_corr_wsum_t$invert_I=ifelse(result_corr_wsum_t$logFC>0,TRUE,FALSE)
        result_corr_wsum_t$P.Value.A=two2one(result_corr_wsum_t$P.Value, two = as.logical(result_corr_wsum_t$two), invert = as.logical(result_corr_wsum_t$invert_A))
        result_corr_wsum_t$P.Value.I=two2one(result_corr_wsum_t$P.Value, two = as.logical(result_corr_wsum_t$two), invert = as.logical(result_corr_wsum_t$invert_I))
        result_corr_wsum_t_A=result_corr_wsum_t[,c("geneset_id","P.Value.A")]
        result_corr_wsum_t_I=result_corr_wsum_t[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wsum_t_A)[2]=paste0(vector[i],".corr_wsum_t")
        colnames(result_corr_wsum_t_I)[2]=paste0(vector[i],".corr_wsum_t")
        rm(result_corr_wsum_t)
      }
      if("corr_wsum_wilcoxon" %in% enrichment_method){
        print("corr_wsum_wilcoxon will start.")
        result=Get_wilcox_test(geneset=geneset,group=group)
        result_corr_wsum_wilcoxon=result[,c("geneset_id","P.Value","logFC")]
        result_corr_wsum_wilcoxon$two=TRUE
        result_corr_wsum_wilcoxon$invert_A=ifelse(result_corr_wsum_wilcoxon$logFC<=0,TRUE,FALSE)
        result_corr_wsum_wilcoxon$invert_I=ifelse(result_corr_wsum_wilcoxon$logFC>0,TRUE,FALSE)
        result_corr_wsum_wilcoxon$P.Value.A=two2one(result_corr_wsum_wilcoxon$P.Value, two = as.logical(result_corr_wsum_wilcoxon$two), invert = as.logical(result_corr_wsum_wilcoxon$invert_A))
        result_corr_wsum_wilcoxon$P.Value.I=two2one(result_corr_wsum_wilcoxon$P.Value, two = as.logical(result_corr_wsum_wilcoxon$two), invert = as.logical(result_corr_wsum_wilcoxon$invert_I))
        result_corr_wsum_wilcoxon_A=result_corr_wsum_wilcoxon[,c("geneset_id","P.Value.A")]
        result_corr_wsum_wilcoxon_I=result_corr_wsum_wilcoxon[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wsum_wilcoxon_A)[2]=paste0(vector[i],".corr_wsum_wilcoxon")
        colnames(result_corr_wsum_wilcoxon_I)[2]=paste0(vector[i],".corr_wsum_wilcoxon")
        rm(result_corr_wsum_wilcoxon)
      }
      if("corr_wsum_kruskal" %in% enrichment_method){
        print("corr_wsum_kruskal will start.")
        result=Get_kruskal_test(geneset=geneset,group=group)
        result_corr_wsum_kruskal=result[,c("geneset_id","P.Value","logFC")]
        result_corr_wsum_kruskal$two=TRUE
        result_corr_wsum_kruskal$invert_A=ifelse(result_corr_wsum_kruskal$logFC<=0,TRUE,FALSE)
        result_corr_wsum_kruskal$invert_I=ifelse(result_corr_wsum_kruskal$logFC>0,TRUE,FALSE)
        result_corr_wsum_kruskal$P.Value.A=two2one(result_corr_wsum_kruskal$P.Value, two = as.logical(result_corr_wsum_kruskal$two), invert = as.logical(result_corr_wsum_kruskal$invert_A))
        result_corr_wsum_kruskal$P.Value.I=two2one(result_corr_wsum_kruskal$P.Value, two = as.logical(result_corr_wsum_kruskal$two), invert = as.logical(result_corr_wsum_kruskal$invert_I))
        result_corr_wsum_kruskal_A=result_corr_wsum_kruskal[,c("geneset_id","P.Value.A")]
        result_corr_wsum_kruskal_I=result_corr_wsum_kruskal[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wsum_kruskal_A)[2]=paste0(vector[i],".corr_wsum_kruskal")
        colnames(result_corr_wsum_kruskal_I)[2]=paste0(vector[i],".corr_wsum_kruskal")
        rm(result_corr_wsum_kruskal)
      }
      if("corr_wsum_permutation" %in% enrichment_method){
        print("corr_wsum_permutation will start.")
        result=Get_permutation_test(geneset=geneset,group=group)
        result_corr_wsum_permutation=result[,c("geneset_id","P.Value","logFC")]
        result_corr_wsum_permutation$two=TRUE
        result_corr_wsum_permutation$invert_A=ifelse(result_corr_wsum_permutation$logFC<=0,TRUE,FALSE)
        result_corr_wsum_permutation$invert_I=ifelse(result_corr_wsum_permutation$logFC>0,TRUE,FALSE)
        result_corr_wsum_permutation=na.omit(result_corr_wsum_permutation)
        result_corr_wsum_permutation$P.Value.A=two2one(result_corr_wsum_permutation$P.Value, two = as.logical(result_corr_wsum_permutation$two), invert = as.logical(result_corr_wsum_permutation$invert_A))
        result_corr_wsum_permutation$P.Value.I=two2one(result_corr_wsum_permutation$P.Value, two = as.logical(result_corr_wsum_permutation$two), invert = as.logical(result_corr_wsum_permutation$invert_I))
        result_corr_wsum_permutation_A=result_corr_wsum_permutation[,c("geneset_id","P.Value.A")]
        result_corr_wsum_permutation_I=result_corr_wsum_permutation[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wsum_permutation_A)[2]=paste0(vector[i],".corr_wsum_permutation")
        colnames(result_corr_wsum_permutation_I)[2]=paste0(vector[i],".corr_wsum_permutation")
        rm(result_corr_wsum_permutation)
      }
      if("corr_wsum_limma" %in% enrichment_method){
        print("corr_wsum_limma will start.")
        result=Get_limma(geneset=geneset,group=group)
        result_corr_wsum_limma=result[,c("geneset_id","P.Value","logFC")]
        result_corr_wsum_limma$two=TRUE
        result_corr_wsum_limma$invert_A=ifelse(result_corr_wsum_limma$logFC<=0,TRUE,FALSE)
        result_corr_wsum_limma$invert_I=ifelse(result_corr_wsum_limma$logFC>0,TRUE,FALSE)
        result_corr_wsum_limma$P.Value.A=two2one(result_corr_wsum_limma$P.Value, two = as.logical(result_corr_wsum_limma$two), invert = as.logical(result_corr_wsum_limma$invert_A))
        result_corr_wsum_limma$P.Value.I=two2one(result_corr_wsum_limma$P.Value, two = as.logical(result_corr_wsum_limma$two), invert = as.logical(result_corr_wsum_limma$invert_I))
        result_corr_wsum_limma_A=result_corr_wsum_limma[,c("geneset_id","P.Value.A")]
        result_corr_wsum_limma_I=result_corr_wsum_limma[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wsum_limma_A)[2]=paste0(vector[i],".corr_wsum_limma")
        colnames(result_corr_wsum_limma_I)[2]=paste0(vector[i],".corr_wsum_limma")
        rm(result_corr_wsum_limma)
      }
      if("corr_wsum_anova" %in% enrichment_method){
        print("corr_wsum_anova will start.")
        result=Get_anova(geneset=geneset,group=group)
        result_corr_wsum_anova=result[,c("geneset_id","P.Value","logFC")]
        result_corr_wsum_anova$two=TRUE
        result_corr_wsum_anova$invert_A=ifelse(result_corr_wsum_anova$logFC<=0,TRUE,FALSE)
        result_corr_wsum_anova$invert_I=ifelse(result_corr_wsum_anova$logFC>0,TRUE,FALSE)
        result_corr_wsum_anova$P.Value.A=two2one(result_corr_wsum_anova$P.Value, two = as.logical(result_corr_wsum_anova$two), invert = as.logical(result_corr_wsum_anova$invert_A))
        result_corr_wsum_anova$P.Value.I=two2one(result_corr_wsum_anova$P.Value, two = as.logical(result_corr_wsum_anova$two), invert = as.logical(result_corr_wsum_anova$invert_I))
        result_corr_wsum_anova_A=result_corr_wsum_anova[,c("geneset_id","P.Value.A")]
        result_corr_wsum_anova_I=result_corr_wsum_anova[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wsum_anova_A)[2]=paste0(vector[i],".corr_wsum_anova")
        colnames(result_corr_wsum_anova_I)[2]=paste0(vector[i],".corr_wsum_anova")
        rm(result_corr_wsum_anova)
      }
    } else {
      print("Not using corr_wsum.")
    }
    merge_result_A=data.frame()
    merge_result_I=data.frame()
    for(m in 1:length(enrichment_method)){
      aggregation_A=get(paste0("result_",enrichment_method[m],"_A"))
      aggregation_I=get(paste0("result_",enrichment_method[m],"_I"))
      aggregation_A$geneset_id=gsub("\\.", " ", aggregation_A$geneset_id)
      aggregation_I$geneset_id=gsub("\\.", " ", aggregation_I$geneset_id)
      if(m==1){
        merge_result_A=aggregation_A
        merge_result_I=aggregation_I
      }else{
        merge_result_A=merge(merge_result_A,aggregation_A,by="geneset_id",all=T)
        merge_result_A[is.na(merge_result_A)]=1
        merge_result_I=merge(merge_result_I,aggregation_I,by="geneset_id",all=T)
        merge_result_I[is.na(merge_result_I)]=1
      }
    }
    assign(paste0(vector[i],"_merge_A"),merge_result_A)
    assign(paste0(vector[i],"_merge_I"),merge_result_I)
  }
  p_matrix_huge_A=data.frame()
  p_matrix_huge_I=data.frame()
  for(datasets in 1:length(vector)){
    p_matrix_method_A=get(paste0(vector[datasets],"_merge_A"))
    p_matrix_method_I=get(paste0(vector[datasets],"_merge_I"))
    if(datasets==1){
      p_matrix_huge_A=p_matrix_method_A
      p_matrix_huge_I=p_matrix_method_I
    }else{
      p_matrix_huge_A=merge(p_matrix_huge_A,p_matrix_method_A,by="geneset_id",all=T)
      p_matrix_huge_A[is.na(p_matrix_huge_A)]=1
      p_matrix_huge_I=merge(p_matrix_huge_I,p_matrix_method_I,by="geneset_id",all=T)
      p_matrix_huge_I[is.na(p_matrix_huge_I)]=1
    }
  }
  rownames(p_matrix_huge_A)=p_matrix_huge_A[,1]
  p_matrix_huge_A=p_matrix_huge_A[,-1,drop=F]
  rownames(p_matrix_huge_I)=p_matrix_huge_I[,1]
  p_matrix_huge_I=p_matrix_huge_I[,-1,drop=F]
  #####################################Calculation of the activation direction index
  final_res_A=data.frame()
  for(inte in 1:length(p_combine_method)){
    col_names=colnames(p_matrix_huge_A)
    prefixes=sapply(strsplit(col_names, "\\."), `[`, 1)
    prefixes_list=split(seq_along(prefixes), prefixes)
    step1=data.frame()
    for (prefix in 1:length(names(prefixes_list))) {
      matrix=prefixes_list[[prefix]]
      matrix=p_matrix_huge_A[,matrix,drop=F]
      if (ncol(matrix) > 1) {
        matrix=as.data.frame(apply(matrix, 2, function(col) {
          if (any(col <= 0, na.rm = TRUE)) {
            min_non_zero=min(col[col > 0], na.rm = TRUE)
            col[col <= 0]=min_non_zero * 0.1
          }
          col[col > 1]=1
          return(col)
        }))
        combined_values=pvalue_integration_strategy(matrix=matrix,p_combine_method=p_combine_method[inte])
        colnames(combined_values)=names(prefixes_list)[prefix]
      }else{
        col=matrix[, 1]
        if (all(col <= 0, na.rm = TRUE)) {
          matrix[, 1]=0
        } else if (any(col <= 0, na.rm = TRUE)) {
          min_non_zero=min(col[col > 0], na.rm = TRUE)
          col[col <= 0]=min_non_zero * 0.1
          matrix[, 1]=col
        }
        matrix[matrix[, 1] > 1, 1]=1
        combined_values=matrix
        colnames(combined_values)=names(prefixes_list)[prefix]
      }
      step1=rbind(step1, as.data.frame(t(combined_values)))
    }
    step1=as.data.frame(t(step1))
    step1=as.data.frame(apply(step1, 2, function(col) {
      if (all(col <= 0, na.rm = TRUE)) {
        col[]=0
      }
      else if (any(col <= 0, na.rm = TRUE)) {
        min_non_zero=min(col[col > 0], na.rm = TRUE)
        col[col <= 0]=min_non_zero * 0.1
      }
      col[col > 1]=1
      return(col)
    }))
    if (all(step1 == 0)) {
      stop("These datasets collected are not suitable for this selected p-value integration method.")
    } else {
      zero_cols = which(colSums(step1) == 0)
      if (length(zero_cols) > 0) {
        colnames_to_remove = colnames(step1)[zero_cols]
        step1 = step1[, -zero_cols, drop = FALSE]
        for (colname in colnames_to_remove) {
          if (grepl("_", colname)) {
            prefix = sub("_.*", "", colname)
            step1 = step1[, !grepl(paste0("^", prefix, "_"), colnames(step1)), drop = FALSE]
          }
        }
      }
      step1 = step1[, colSums(step1) > 0, drop = FALSE]
    }
    col_names=colnames(step1)
    prefixes = sapply(strsplit(col_names, "_"), `[`, 1)
    prefixes_list=split(seq_along(prefixes), prefixes)
    step2=data.frame()
    for (prefix in 1:length(names(prefixes_list))) {
      matrix=prefixes_list[[prefix]]
      matrix=step1[,matrix,drop=F]
      if (ncol(matrix) > 1) {
        matrix=as.data.frame(apply(matrix, 2, function(col) {
          if (any(col <= 0, na.rm = TRUE)) {
            min_non_zero=min(col[col > 0], na.rm = TRUE)
            col[col <= 0]=min_non_zero * 0.1
          }
          col[col > 1]=1
          return(col)
        }))
        combined_values=pvalue_integration_strategy(matrix=matrix,p_combine_method=p_combine_method[inte])
        colnames(combined_values)=names(prefixes_list)[prefix]
      }else{
        col=matrix[, 1]
        if (all(col <= 0, na.rm = TRUE)) {
          matrix[, 1]=0
        } else if (any(col <= 0, na.rm = TRUE)) {
          min_non_zero=min(col[col > 0], na.rm = TRUE)
          col[col <= 0]=min_non_zero * 0.1
          matrix[, 1]=col
        }
        matrix[matrix[, 1] > 1, 1]=1
        combined_values=matrix
        colnames(combined_values)=names(prefixes_list)[prefix]
      }
      step2=rbind(step2, as.data.frame(t(combined_values)))
    }
    step2=as.data.frame(t(step2))
    step2 = as.data.frame(apply(step2, 2, function(col) {
      if (all(col <= 0, na.rm = TRUE)) {
        col[] = 0
      } else if (any(col <= 0, na.rm = TRUE)) {
        min_non_zero = min(col[col > 0], na.rm = TRUE)
        col[col <= 0] = min_non_zero * 0.1
      }
      col[col > 1]=1
      return(col)
    }))
    if (all(step2 == 0)) {
      stop("These datasets collected are not suitable for this selected p-value integration method.")
    } else {
      step2 = step2[, colSums(step2) > 0, drop = FALSE]
    }
    step2 = as.data.frame(apply(step2, 2, function(col) {
      min_non_zero = min(col[col > 0], na.rm = TRUE)
      if (!is.infinite(min_non_zero)) {
        col[col <= 0] = min_non_zero
      }
      return(col)
    }))
    step2_value=as.data.frame(rowMeans(-log10(step2)))
    colnames(step2_value)=p_combine_method[inte]
    step2_value=as.data.frame(t(step2_value))
    final_res_A=rbind(final_res_A,step2_value)
  }
  final_res_A=as.data.frame(t(final_res_A))
  final_res_A=as.data.frame(rowMeans(final_res_A))
  colnames(final_res_A)="ADI"
  #####################################Calculation of the inhibition direction index
  final_res_I=data.frame()
  for(inte in 1:length(p_combine_method)){
    col_names=colnames(p_matrix_huge_I)
    prefixes=sapply(strsplit(col_names, "\\."), `[`, 1)
    prefixes_list=split(seq_along(prefixes), prefixes)
    step1=data.frame()
    for (prefix in 1:length(names(prefixes_list))) {
      matrix=prefixes_list[[prefix]]
      matrix=p_matrix_huge_I[,matrix,drop=F]
      if (ncol(matrix) > 1) {
        matrix=as.data.frame(apply(matrix, 2, function(col) {
          if (any(col <= 0, na.rm = TRUE)) {
            min_non_zero=min(col[col > 0], na.rm = TRUE)
            col[col <= 0]=min_non_zero * 0.1
          }
          col[col > 1]=1
          return(col)
        }))
        combined_values=pvalue_integration_strategy(matrix=matrix,p_combine_method=p_combine_method[inte])
        colnames(combined_values)=names(prefixes_list)[prefix]
      }else{
        col=matrix[, 1]
        if (all(col <= 0, na.rm = TRUE)) {
          matrix[, 1]=0
        } else if (any(col <= 0, na.rm = TRUE)) {
          min_non_zero=min(col[col > 0], na.rm = TRUE)
          col[col <= 0]=min_non_zero * 0.1
          matrix[, 1]=col
        }
        matrix[matrix[, 1] > 1, 1]=1
        combined_values=matrix
        colnames(combined_values)=names(prefixes_list)[prefix]
      }
      step1=rbind(step1, as.data.frame(t(combined_values)))
    }
    step1=as.data.frame(t(step1))
    step1=as.data.frame(apply(step1, 2, function(col) {
      if (all(col <= 0, na.rm = TRUE)) {
        col[]=0
      }
      else if (any(col <= 0, na.rm = TRUE)) {
        min_non_zero=min(col[col > 0], na.rm = TRUE)
        col[col <= 0]=min_non_zero * 0.1
      }
      col[col > 1]=1
      return(col)
    }))
    if (all(step1 == 0)) {
      stop("These datasets collected are not suitable for this selected p-value integration method.")
    } else {
      zero_cols = which(colSums(step1) == 0)
      if (length(zero_cols) > 0) {
        colnames_to_remove = colnames(step1)[zero_cols]
        step1 = step1[, -zero_cols, drop = FALSE]
        for (colname in colnames_to_remove) {
          if (grepl("_", colname)) {
            prefix = sub("_.*", "", colname)
            step1 = step1[, !grepl(paste0("^", prefix, "_"), colnames(step1)), drop = FALSE]
          }
        }
      }
      step1 = step1[, colSums(step1) > 0, drop = FALSE]
    }
    col_names=colnames(step1)
    prefixes = sapply(strsplit(col_names, "_"), `[`, 1)
    prefixes_list=split(seq_along(prefixes), prefixes)
    step2=data.frame()
    for (prefix in 1:length(names(prefixes_list))) {
      matrix=prefixes_list[[prefix]]
      matrix=step1[,matrix,drop=F]
      if (ncol(matrix) > 1) {
        matrix=as.data.frame(apply(matrix, 2, function(col) {
          if (any(col <= 0, na.rm = TRUE)) {
            min_non_zero=min(col[col > 0], na.rm = TRUE)
            col[col <= 0]=min_non_zero * 0.1
          }
          col[col > 1]=1
          return(col)
        }))
        combined_values=pvalue_integration_strategy(matrix=matrix,p_combine_method=p_combine_method[inte])
        colnames(combined_values)=names(prefixes_list)[prefix]
      }else{
        col=matrix[, 1]
        if (all(col <= 0, na.rm = TRUE)) {
          matrix[, 1]=0
        } else if (any(col <= 0, na.rm = TRUE)) {
          min_non_zero=min(col[col > 0], na.rm = TRUE)
          col[col <= 0]=min_non_zero * 0.1
          matrix[, 1]=col
        }
        matrix[matrix[, 1] > 1, 1]=1
        combined_values=matrix
        colnames(combined_values)=names(prefixes_list)[prefix]
      }
      step2=rbind(step2, as.data.frame(t(combined_values)))
    }
    step2=as.data.frame(t(step2))
    step2 = as.data.frame(apply(step2, 2, function(col) {
      if (all(col <= 0, na.rm = TRUE)) {
        col[] = 0
      } else if (any(col <= 0, na.rm = TRUE)) {
        min_non_zero = min(col[col > 0], na.rm = TRUE)
        col[col <= 0] = min_non_zero * 0.1
      }
      col[col > 1]=1
      return(col)
    }))
    if (all(step2 == 0)) {
      stop("These datasets collected are not suitable for this selected p-value integration method.")
    } else {
      step2 = step2[, colSums(step2) > 0, drop = FALSE]
    }
    step2 = as.data.frame(apply(step2, 2, function(col) {
      min_non_zero = min(col[col > 0], na.rm = TRUE)
      if (!is.infinite(min_non_zero)) {
        col[col <= 0] = min_non_zero
      }
      return(col)
    }))
    step2_value=as.data.frame(rowMeans(log10(step2)))
    colnames(step2_value)=p_combine_method[inte]
    step2_value=as.data.frame(t(step2_value))
    final_res_I=rbind(final_res_I,step2_value)
  }
  final_res_I=as.data.frame(t(final_res_I))
  final_res_I=as.data.frame(rowMeans(final_res_I))
  colnames(final_res_I)="IDI"
  final_res=merge(final_res_A,final_res_I,by="row.names",all=T)
  final_res[is.na(final_res)]=0
  colnames(final_res)[1]="signatures"
  if (!setequal(label_genesets, final_res$signatures)) {
    new_signatures=setdiff(label_genesets, final_res$signatures)
    if (length(new_signatures) > 0) {
      new_rows=data.frame(signatures = new_signatures, ADI = rep(0, length(new_signatures)), IDI = rep(0, length(new_signatures)))
      final_res=rbind(final_res, new_rows)  
    }
  }
  return(final_res)
}