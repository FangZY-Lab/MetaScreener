#' @title CorMetaScreener: a multi-level meta-analysis framework robustly screens directional signatures based on correlation analyses.
#' @param expression_accession_vector A vector of all the sub-dataset names.
#' @param geneSets_gmt Genesets to be screened.(Please make sure that every gene set in there is unique, otherwise there will be severe covariance.)
#' @param enrichment_method Selection of enrichment methods, multiple options available.
#' @param min.sz The genesets contain the minimum number of genes, where the default parameter is 1.
#' @param max.sz The genesets contain the maximum number of genes, where the default parameter is 10000.
#' @param p_combine_method The choice of method for integrating p-values.
#' @author Dingkang Zhao
#' @examples
#' @return
#' @export
CorMetaScreener=function(expression_accession_vector,
                          geneSets_gmt,
                          enrichment_method=c("gsva_pearson","gsva_kendall","gsva_spearman","gsva_lm","gsva_biweight","gsva_distance","gsva_percentage","gsva_blomqvist","gsva_hoeffding","gsva_gamma",
                                              "ssgsea_pearson","ssgsea_kendall","ssgsea_spearman","ssgsea_lm","ssgsea_biweight","ssgsea_distance","ssgsea_percentage","ssgsea_blomqvist","ssgsea_hoeffding","ssgsea_gamma",
                                              "zscore_pearson","zscore_kendall","zscore_spearman","zscore_lm","zscore_biweight","zscore_distance","zscore_percentage","zscore_blomqvist","zscore_hoeffding","zscore_gamma",
                                              "plage_pearson","plage_kendall","plage_spearman","plage_lm","plage_biweight","plage_distance","plage_percentage","plage_blomqvist","plage_hoeffding","plage_gamma",
                                              "pca_pearson","pca_kendall","pca_spearman","pca_lm","pca_biweight","pca_distance","pca_percentage","pca_blomqvist","pca_hoeffding","pca_gamma",
                                              "aucell_pearson","aucell_kendall","aucell_spearman","aucell_lm","aucell_biweight","aucell_distance","aucell_percentage","aucell_blomqvist","aucell_hoeffding","aucell_gamma",
                                              "ucell_pearson","ucell_kendall","ucell_spearman","ucell_lm","ucell_biweight","ucell_distance","ucell_percentage","ucell_blomqvist","ucell_hoeffding","ucell_gamma",
                                              "singscore_pearson","singscore_kendall","singscore_spearman","singscore_lm","singscore_biweight","singscore_distance","singscore_percentage","singscore_blomqvist","singscore_hoeffding","singscore_gamma",
                                              "median_pearson","median_kendall","median_spearman","median_lm","median_biweight","median_distance","median_percentage","median_blomqvist","median_hoeffding","median_gamma",
                                              "pearson_fgsea","kendall_fgsea","spearman_fgsea","lm_fgsea","biweight_fgsea","distance_fgsea","percentage_fgsea","blomqvist_fgsea","hoeffding_fgsea","gamma_fgsea",
                                              "pearson_ora","kendall_ora","spearman_ora","lm_ora","biweight_ora","distance_ora","percentage_ora","blomqvist_ora","hoeffding_ora","gamma_ora",
                                              "consensus_pearson","consensus_kendall","consensus_spearman","consensus_lm","consensus_biweight","consensus_distance","consensus_percentage","consensus_blomqvist","consensus_hoeffding","consensus_gamma",
                                              "mdt_pearson","mdt_kendall","mdt_spearman","mdt_lm","mdt_biweight","mdt_distance","mdt_percentage","mdt_blomqvist","mdt_hoeffding","mdt_gamma",
                                              "mlm_pearson","mlm_kendall","mlm_spearman","mlm_lm","mlm_biweight","mlm_distance","mlm_percentage","mlm_blomqvist","mlm_hoeffding","mlm_gamma",
                                              "udt_pearson","udt_kendall","udt_spearman","udt_lm","udt_biweight","udt_distance","udt_percentage","udt_blomqvist","udt_hoeffding","udt_gamma",
                                              "ulm_pearson","ulm_kendall","ulm_spearman","ulm_lm","ulm_biweight","ulm_distance","ulm_percentage","ulm_blomqvist","ulm_hoeffding","ulm_gamma",
                                              "viper_pearson","viper_kendall","viper_spearman","viper_lm","viper_biweight","viper_distance","viper_percentage","viper_blomqvist","viper_hoeffding","viper_gamma",
                                              "wmean_pearson","wmean_kendall","wmean_spearman","wmean_lm","wmean_biweight","wmean_distance","wmean_percentage","wmean_blomqvist","wmean_hoeffding","wmean_gamma",
                                              "norm_wmean_pearson","norm_wmean_kendall","norm_wmean_spearman","norm_wmean_lm","norm_wmean_biweight","norm_wmean_distance","norm_wmean_percentage","norm_wmean_blomqvist","norm_wmean_hoeffding","norm_wmean_gamma",
                                              "corr_wmean_pearson","corr_wmean_kendall","corr_wmean_spearman","corr_wmean_lm","corr_wmean_biweight","corr_wmean_distance","corr_wmean_percentage","corr_wmean_blomqvist","corr_wmean_hoeffding","corr_wmean_gamma",
                                              "wsum_pearson","wsum_kendall","wsum_spearman","wsum_lm","wsum_biweight","wsum_distance","wsum_percentage","wsum_blomqvist","wsum_hoeffding","wsum_gamma",
                                              "norm_wsum_pearson","norm_wsum_kendall","norm_wsum_spearman","norm_wsum_lm","norm_wsum_biweight","norm_wsum_distance","norm_wsum_percentage","norm_wsum_blomqvist","norm_wsum_hoeffding","norm_wsum_gamma",
                                              "corr_wsum_pearson","corr_wsum_kendall","corr_wsum_spearman","corr_wsum_lm","corr_wsum_biweight","corr_wsum_distance","corr_wsum_percentage","corr_wsum_blomqvist","corr_wsum_hoeffding","corr_wsum_gamma"),
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
  library(plyr)
  library(reshape2)
  library(GSEABase)
  library(GSVA)
  library(doBy)
  library(dplyr)
  library(tidyverse)
  library(devtools)
  library(stringr)
  library(correlation)
  library(wdm)
  library(caret)
  library(fgsea)
  library(viper)
  ##########################################CTT: a p-value integration method.
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
  ##########################################pearson_test
  pearson_test=function(geneset,value){
    geneset=as.data.frame(t(geneset))
    id=colnames(geneset)
    res_all=data.frame()
    for(edu in 1:length(id)){
      R=cor.test(value$value,geneset[,edu],method = "pearson",alternative="two.sided")
      p.value=R$p.value
      relation=R$estimate
      res=data.frame(geneset_id=id[edu],P.Value=p.value,r=relation)
      res_all=rbind(res_all,res)
    }
    row.names(res_all)=NULL
    return(res_all)
  }
  ##########################################kendall_test
  kendall_test=function(geneset,value){
    geneset=as.data.frame(t(geneset))
    id=colnames(geneset)
    res_all=data.frame()
    for(edu in 1:length(id)){
      R=cor.test(value$value,geneset[,edu],method = "kendall",alternative="two.sided")
      p.value=R$p.value
      relation=R$estimate
      res=data.frame(geneset_id=id[edu],P.Value=p.value,r=relation)
      res_all=rbind(res_all,res)
    }
    row.names(res_all)=NULL
    return(res_all)
  }
  ##########################################spearman_test
  spearman_test=function(geneset,value){
    geneset=as.data.frame(t(geneset))
    id=colnames(geneset)
    res_all=data.frame()
    for(edu in 1:length(id)){
      R=cor.test(value$value,geneset[,edu],method = "spearman",alternative="two.sided")
      p.value=R$p.value
      relation=R$estimate
      res=data.frame(geneset_id=id[edu],P.Value=p.value,r=relation)
      res_all=rbind(res_all,res)
    }
    row.names(res_all)=NULL
    return(res_all)
  }
  ##########################################lm_test
  lm_test=function(geneset,value){
    geneset=as.data.frame(t(geneset))
    id=colnames(geneset)
    res_all=data.frame()
    for(edu in 1:length(id)){
      R=lm(value$value~geneset[,edu])
      R=summary(R)
      if (nrow(R$coefficients) >= 2) {
        p.value = R$coefficients[2, 4]
        relation = R$coefficients[2, 1]
      } else {
        p.value=1
        relation=0
      }
      res=data.frame(geneset_id=id[edu],P.Value=p.value,r=relation)
      res_all=rbind(res_all,res)
    }
    row.names(res_all)=NULL
    return(res_all)
  }
  ##########################################biweight_test
  biweight_test=function(geneset,value){
    geneset_mm=as.data.frame(t(geneset))
    res_all=data.frame()
    for(mm in 1:length(colnames(geneset_mm))){
      res=correlation(value$value,geneset_mm[,mm,drop=F],p_adjust="none",method="biweight")
      res=as.data.frame(res[,c(2,9,3)])
      res_all=rbind(res_all,res)
    }
    colnames(res_all)=c("geneset_id","P.Value","r")
    res_all$P.Value[is.na(res_all$P.Value)]=1
    res_all$r[is.na(res_all$r)]=0
    return(res_all)
  }
  ##########################################distance_test
  distance_test=function(geneset,value){
    geneset_mm=as.data.frame(t(geneset))
    res_all=data.frame()
    for(mm in 1:length(colnames(geneset_mm))){
      res=correlation(value$value,geneset_mm[,mm,drop=F],p_adjust="none",method="distance")
      res=as.data.frame(res[,c(2,9,3)])
      res_all=rbind(res_all,res)
    }
    colnames(res_all)=c("geneset_id","P.Value","r")
    res_all$P.Value[is.na(res_all$P.Value)]=1
    res_all$r[is.na(res_all$r)]=0
    return(res_all)
  }
  ##########################################percentage_test
  percentage_test=function(geneset,value){
    geneset_mm=as.data.frame(t(geneset))
    res_all=data.frame()
    for(mm in 1:length(colnames(geneset_mm))){
      res=correlation(value$value,geneset_mm[,mm,drop=F],p_adjust="none",method="percentage")
      res=as.data.frame(res[,c(2,9,3)])
      res_all=rbind(res_all,res)
    }
    colnames(res_all)=c("geneset_id","P.Value","r")
    res_all$P.Value[is.na(res_all$P.Value)]=1
    res_all$r[is.na(res_all$r)]=0
    return(res_all)
  }
  ##########################################blomqvist_test
  blomqvist_test=function(geneset,value){
    geneset_mm=as.data.frame(t(geneset))
    res_all=data.frame()
    for(mm in 1:length(colnames(geneset_mm))){
      res=correlation(value$value,geneset_mm[,mm,drop=F],p_adjust="none",method="blomqvist")
      res=as.data.frame(res[,c(2,9,3)])
      res_all=rbind(res_all,res)
    }
    colnames(res_all)=c("geneset_id","P.Value","r")
    res_all$P.Value[is.na(res_all$P.Value)]=1
    res_all$r[is.na(res_all$r)]=0
    return(res_all)
  }
  ##########################################hoeffding_test
  hoeffding_test=function(geneset,value){
    geneset_mm=as.data.frame(t(geneset))
    res_all=data.frame()
    for(mm in 1:length(colnames(geneset_mm))){
      res=correlation(value$value,geneset_mm[,mm,drop=F],p_adjust="none",method="hoeffding")
      res=as.data.frame(res[,c(2,9,3)])
      res_all=rbind(res_all,res)
    }
    colnames(res_all)=c("geneset_id","P.Value","r")
    res_all$P.Value[is.na(res_all$P.Value)]=1
    res_all$r[is.na(res_all$r)]=0
    return(res_all)
  }
  ##########################################gamma_test
  gamma_test=function(geneset,value){
    geneset_mm=as.data.frame(t(geneset))
    res_all=data.frame()
    for(mm in 1:length(colnames(geneset_mm))){
      res=correlation(value$value,geneset_mm[,mm,drop=F],p_adjust="none",method="gamma")
      res=as.data.frame(res[,c(2,9,3)])
      res_all=rbind(res_all,res)
    }
    colnames(res_all)=c("geneset_id","P.Value","r")
    res_all$P.Value[is.na(res_all$P.Value)]=1
    res_all$r[is.na(res_all$r)]=0
    return(res_all)
  }
  ##########################################
  vector=expression_accession_vector
  print("Start using multiple enrichment analysis methods to cycle through the calculation of enrichment scores for each gene set in each dataset.")
  for (i in 1:length(vector)){
    print("Start loading expression data and numeric information.")
    print(paste0(vector[i]))
    exp=get(vector[i])
    dimnames=list(rownames(exp),colnames(exp))
    exp=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
    exp=as.data.frame(exp)
    constant_cols=sapply(exp, function(x) var(x) == 0)
    exp=exp[, !constant_cols]
    constant_rows=apply(exp, 1, function(x) length(unique(x)) == 1)
    exp=exp[!constant_rows, ]
    value=get(paste0(vector[i],"_V"))
    value=as.data.frame(value)
    geneSets=geneSets_gmt
    label_genesets=names(geneSets)
    ##########################################################################Selection of the method of calculation
    ##########################################gsva
    methods_to_check=c("gsva_pearson","gsva_kendall","gsva_spearman","gsva_lm","gsva_biweight","gsva_distance","gsva_percentage","gsva_blomqvist","gsva_hoeffding","gsva_gamma")
    if (any(methods_to_check %in% enrichment_method)){
      geneset=gsvaParam(exprData=as.matrix(exp),geneSets=geneSets,kcdf="Gaussian",minSize=min.sz,maxSize=max.sz)
      geneset=gsva(geneset)
      geneset=as.data.frame(geneset)
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("gsva_pearson" %in% enrichment_method){
        print("gsva_pearson will start.")
        result_gsva_pearson=pearson_test(geneset=geneset,value=value)
        result_gsva_pearson$two=TRUE
        result_gsva_pearson$invert_A=ifelse(result_gsva_pearson$r<=0,TRUE,FALSE)
        result_gsva_pearson$invert_I=ifelse(result_gsva_pearson$r>0,TRUE,FALSE)
        result_gsva_pearson$P.Value.A=two2one(result_gsva_pearson$P.Value, two = as.logical(result_gsva_pearson$two), invert = as.logical(result_gsva_pearson$invert_A))
        result_gsva_pearson$P.Value.I=two2one(result_gsva_pearson$P.Value, two = as.logical(result_gsva_pearson$two), invert = as.logical(result_gsva_pearson$invert_I))
        result_gsva_pearson_A=result_gsva_pearson[,c("geneset_id","P.Value.A")]
        result_gsva_pearson_I=result_gsva_pearson[,c("geneset_id","P.Value.I")]
        colnames(result_gsva_pearson_A)[2]=paste0(vector[i],".gsva_pearson")
        colnames(result_gsva_pearson_I)[2]=paste0(vector[i],".gsva_pearson")
        rm(result_gsva_pearson)
      }
      if("gsva_kendall" %in% enrichment_method){
        print("gsva_kendall will start.")
        result_gsva_kendall=kendall_test(geneset=geneset,value=value)
        result_gsva_kendall$two=TRUE
        result_gsva_kendall$invert_A=ifelse(result_gsva_kendall$r<=0,TRUE,FALSE)
        result_gsva_kendall$invert_I=ifelse(result_gsva_kendall$r>0,TRUE,FALSE)
        result_gsva_kendall$P.Value.A=two2one(result_gsva_kendall$P.Value, two = as.logical(result_gsva_kendall$two), invert = as.logical(result_gsva_kendall$invert_A))
        result_gsva_kendall$P.Value.I=two2one(result_gsva_kendall$P.Value, two = as.logical(result_gsva_kendall$two), invert = as.logical(result_gsva_kendall$invert_I))
        result_gsva_kendall_A=result_gsva_kendall[,c("geneset_id","P.Value.A")]
        result_gsva_kendall_I=result_gsva_kendall[,c("geneset_id","P.Value.I")]
        colnames(result_gsva_kendall_A)[2]=paste0(vector[i],".gsva_kendall")
        colnames(result_gsva_kendall_I)[2]=paste0(vector[i],".gsva_kendall")
        rm(result_gsva_kendall)
      }
      if("gsva_spearman" %in% enrichment_method){
        print("gsva_spearman will start.")
        result_gsva_spearman=spearman_test(geneset=geneset,value=value)
        result_gsva_spearman$two=TRUE
        result_gsva_spearman$invert_A=ifelse(result_gsva_spearman$r<=0,TRUE,FALSE)
        result_gsva_spearman$invert_I=ifelse(result_gsva_spearman$r>0,TRUE,FALSE)
        result_gsva_spearman$P.Value.A=two2one(result_gsva_spearman$P.Value, two = as.logical(result_gsva_spearman$two), invert = as.logical(result_gsva_spearman$invert_A))
        result_gsva_spearman$P.Value.I=two2one(result_gsva_spearman$P.Value, two = as.logical(result_gsva_spearman$two), invert = as.logical(result_gsva_spearman$invert_I))
        result_gsva_spearman_A=result_gsva_spearman[,c("geneset_id","P.Value.A")]
        result_gsva_spearman_I=result_gsva_spearman[,c("geneset_id","P.Value.I")]
        colnames(result_gsva_spearman_A)[2]=paste0(vector[i],".gsva_spearman")
        colnames(result_gsva_spearman_I)[2]=paste0(vector[i],".gsva_spearman")
        rm(result_gsva_spearman)
      }
      if("gsva_lm" %in% enrichment_method){
        print("gsva_lm will start.")
        result_gsva_lm=lm_test(geneset=geneset,value=value)
        result_gsva_lm$two=TRUE
        result_gsva_lm$invert_A=ifelse(result_gsva_lm$r<=0,TRUE,FALSE)
        result_gsva_lm$invert_I=ifelse(result_gsva_lm$r>0,TRUE,FALSE)
        result_gsva_lm$P.Value.A=two2one(result_gsva_lm$P.Value, two = as.logical(result_gsva_lm$two), invert = as.logical(result_gsva_lm$invert_A))
        result_gsva_lm$P.Value.I=two2one(result_gsva_lm$P.Value, two = as.logical(result_gsva_lm$two), invert = as.logical(result_gsva_lm$invert_I))
        result_gsva_lm_A=result_gsva_lm[,c("geneset_id","P.Value.A")]
        result_gsva_lm_I=result_gsva_lm[,c("geneset_id","P.Value.I")]
        colnames(result_gsva_lm_A)[2]=paste0(vector[i],".gsva_lm")
        colnames(result_gsva_lm_I)[2]=paste0(vector[i],".gsva_lm")
        rm(result_gsva_lm)
      }
      if("gsva_biweight" %in% enrichment_method){
        print("gsva_biweight will start.")
        result_gsva_biweight=biweight_test(geneset=geneset,value=value)
        result_gsva_biweight$two=TRUE
        result_gsva_biweight$invert_A=ifelse(result_gsva_biweight$r<=0,TRUE,FALSE)
        result_gsva_biweight$invert_I=ifelse(result_gsva_biweight$r>0,TRUE,FALSE)
        result_gsva_biweight$P.Value.A=two2one(result_gsva_biweight$P.Value, two = as.logical(result_gsva_biweight$two), invert = as.logical(result_gsva_biweight$invert_A))
        result_gsva_biweight$P.Value.I=two2one(result_gsva_biweight$P.Value, two = as.logical(result_gsva_biweight$two), invert = as.logical(result_gsva_biweight$invert_I))
        result_gsva_biweight_A=result_gsva_biweight[,c("geneset_id","P.Value.A")]
        result_gsva_biweight_I=result_gsva_biweight[,c("geneset_id","P.Value.I")]
        colnames(result_gsva_biweight_A)[2]=paste0(vector[i],".gsva_biweight")
        colnames(result_gsva_biweight_I)[2]=paste0(vector[i],".gsva_biweight")
        rm(result_gsva_biweight)
      }
      if("gsva_distance" %in% enrichment_method){
        print("gsva_distance will start.")
        result_gsva_distance=distance_test(geneset=geneset,value=value)
        result_gsva_distance$two=TRUE
        result_gsva_distance$invert_A=ifelse(result_gsva_distance$r<=0,TRUE,FALSE)
        result_gsva_distance$invert_I=ifelse(result_gsva_distance$r>0,TRUE,FALSE)
        result_gsva_distance$P.Value.A=two2one(result_gsva_distance$P.Value, two = as.logical(result_gsva_distance$two), invert = as.logical(result_gsva_distance$invert_A))
        result_gsva_distance$P.Value.I=two2one(result_gsva_distance$P.Value, two = as.logical(result_gsva_distance$two), invert = as.logical(result_gsva_distance$invert_I))
        result_gsva_distance_A=result_gsva_distance[,c("geneset_id","P.Value.A")]
        result_gsva_distance_I=result_gsva_distance[,c("geneset_id","P.Value.I")]
        colnames(result_gsva_distance_A)[2]=paste0(vector[i],".gsva_distance")
        colnames(result_gsva_distance_I)[2]=paste0(vector[i],".gsva_distance")
        rm(result_gsva_distance)
      }
      if("gsva_percentage" %in% enrichment_method){
        print("gsva_percentage will start.")
        result_gsva_percentage=percentage_test(geneset=geneset,value=value)
        result_gsva_percentage$two=TRUE
        result_gsva_percentage$invert_A=ifelse(result_gsva_percentage$r<=0,TRUE,FALSE)
        result_gsva_percentage$invert_I=ifelse(result_gsva_percentage$r>0,TRUE,FALSE)
        result_gsva_percentage$P.Value.A=two2one(result_gsva_percentage$P.Value, two = as.logical(result_gsva_percentage$two), invert = as.logical(result_gsva_percentage$invert_A))
        result_gsva_percentage$P.Value.I=two2one(result_gsva_percentage$P.Value, two = as.logical(result_gsva_percentage$two), invert = as.logical(result_gsva_percentage$invert_I))
        result_gsva_percentage_A=result_gsva_percentage[,c("geneset_id","P.Value.A")]
        result_gsva_percentage_I=result_gsva_percentage[,c("geneset_id","P.Value.I")]
        colnames(result_gsva_percentage_A)[2]=paste0(vector[i],".gsva_percentage")
        colnames(result_gsva_percentage_I)[2]=paste0(vector[i],".gsva_percentage")
        rm(result_gsva_percentage)
      }
      if("gsva_blomqvist" %in% enrichment_method){
        print("gsva_blomqvist will start.")
        result_gsva_blomqvist=blomqvist_test(geneset=geneset,value=value)
        result_gsva_blomqvist$two=TRUE
        result_gsva_blomqvist$invert_A=ifelse(result_gsva_blomqvist$r<=0,TRUE,FALSE)
        result_gsva_blomqvist$invert_I=ifelse(result_gsva_blomqvist$r>0,TRUE,FALSE)
        result_gsva_blomqvist$P.Value.A=two2one(result_gsva_blomqvist$P.Value, two = as.logical(result_gsva_blomqvist$two), invert = as.logical(result_gsva_blomqvist$invert_A))
        result_gsva_blomqvist$P.Value.I=two2one(result_gsva_blomqvist$P.Value, two = as.logical(result_gsva_blomqvist$two), invert = as.logical(result_gsva_blomqvist$invert_I))
        result_gsva_blomqvist_A=result_gsva_blomqvist[,c("geneset_id","P.Value.A")]
        result_gsva_blomqvist_I=result_gsva_blomqvist[,c("geneset_id","P.Value.I")]
        colnames(result_gsva_blomqvist_A)[2]=paste0(vector[i],".gsva_blomqvist")
        colnames(result_gsva_blomqvist_I)[2]=paste0(vector[i],".gsva_blomqvist")
        rm(result_gsva_blomqvist)
      }
      if("gsva_hoeffding" %in% enrichment_method){
        print("gsva_hoeffding will start.")
        result_gsva_hoeffding=hoeffding_test(geneset=geneset,value=value)
        result_gsva_hoeffding$two=TRUE
        result_gsva_hoeffding$invert_A=ifelse(result_gsva_hoeffding$r<=0,TRUE,FALSE)
        result_gsva_hoeffding$invert_I=ifelse(result_gsva_hoeffding$r>0,TRUE,FALSE)
        result_gsva_hoeffding$P.Value.A=two2one(result_gsva_hoeffding$P.Value, two = as.logical(result_gsva_hoeffding$two), invert = as.logical(result_gsva_hoeffding$invert_A))
        result_gsva_hoeffding$P.Value.I=two2one(result_gsva_hoeffding$P.Value, two = as.logical(result_gsva_hoeffding$two), invert = as.logical(result_gsva_hoeffding$invert_I))
        result_gsva_hoeffding_A=result_gsva_hoeffding[,c("geneset_id","P.Value.A")]
        result_gsva_hoeffding_I=result_gsva_hoeffding[,c("geneset_id","P.Value.I")]
        colnames(result_gsva_hoeffding_A)[2]=paste0(vector[i],".gsva_hoeffding")
        colnames(result_gsva_hoeffding_I)[2]=paste0(vector[i],".gsva_hoeffding")
        rm(result_gsva_hoeffding)
      }
      if("gsva_gamma" %in% enrichment_method){
        print("gsva_gamma will start.")
        result_gsva_gamma=gamma_test(geneset=geneset,value=value)
        result_gsva_gamma$two=TRUE
        result_gsva_gamma$invert_A=ifelse(result_gsva_gamma$r<=0,TRUE,FALSE)
        result_gsva_gamma$invert_I=ifelse(result_gsva_gamma$r>0,TRUE,FALSE)
        result_gsva_gamma$P.Value.A=two2one(result_gsva_gamma$P.Value, two = as.logical(result_gsva_gamma$two), invert = as.logical(result_gsva_gamma$invert_A))
        result_gsva_gamma$P.Value.I=two2one(result_gsva_gamma$P.Value, two = as.logical(result_gsva_gamma$two), invert = as.logical(result_gsva_gamma$invert_I))
        result_gsva_gamma_A=result_gsva_gamma[,c("geneset_id","P.Value.A")]
        result_gsva_gamma_I=result_gsva_gamma[,c("geneset_id","P.Value.I")]
        colnames(result_gsva_gamma_A)[2]=paste0(vector[i],".gsva_gamma")
        colnames(result_gsva_gamma_I)[2]=paste0(vector[i],".gsva_gamma")
        rm(result_gsva_gamma)
      }
    } else {
      print("Not using gsva.")
    }
    ##########################################ssgsea
    methods_to_check=c("ssgsea_pearson","ssgsea_kendall","ssgsea_spearman","ssgsea_lm","ssgsea_biweight","ssgsea_distance","ssgsea_percentage","ssgsea_blomqvist","ssgsea_hoeffding","ssgsea_gamma")
    if (any(methods_to_check %in% enrichment_method)){
      geneset=ssgseaParam(exprData=as.matrix(exp),geneSets=geneSets,minSize=min.sz,maxSize=max.sz)
      geneset=gsva(geneset)
      geneset=as.data.frame(geneset)
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("ssgsea_pearson" %in% enrichment_method){
        print("ssgsea_pearson will start.")
        result_ssgsea_pearson=pearson_test(geneset=geneset,value=value)
        result_ssgsea_pearson$two=TRUE
        result_ssgsea_pearson$invert_A=ifelse(result_ssgsea_pearson$r<=0,TRUE,FALSE)
        result_ssgsea_pearson$invert_I=ifelse(result_ssgsea_pearson$r>0,TRUE,FALSE)
        result_ssgsea_pearson$P.Value.A=two2one(result_ssgsea_pearson$P.Value, two = as.logical(result_ssgsea_pearson$two), invert = as.logical(result_ssgsea_pearson$invert_A))
        result_ssgsea_pearson$P.Value.I=two2one(result_ssgsea_pearson$P.Value, two = as.logical(result_ssgsea_pearson$two), invert = as.logical(result_ssgsea_pearson$invert_I))
        result_ssgsea_pearson_A=result_ssgsea_pearson[,c("geneset_id","P.Value.A")]
        result_ssgsea_pearson_I=result_ssgsea_pearson[,c("geneset_id","P.Value.I")]
        colnames(result_ssgsea_pearson_A)[2]=paste0(vector[i],".ssgsea_pearson")
        colnames(result_ssgsea_pearson_I)[2]=paste0(vector[i],".ssgsea_pearson")
        rm(result_ssgsea_pearson)
      }
      if("ssgsea_kendall" %in% enrichment_method){
        print("ssgsea_kendall will start.")
        result_ssgsea_kendall=kendall_test(geneset=geneset,value=value)
        result_ssgsea_kendall$two=TRUE
        result_ssgsea_kendall$invert_A=ifelse(result_ssgsea_kendall$r<=0,TRUE,FALSE)
        result_ssgsea_kendall$invert_I=ifelse(result_ssgsea_kendall$r>0,TRUE,FALSE)
        result_ssgsea_kendall$P.Value.A=two2one(result_ssgsea_kendall$P.Value, two = as.logical(result_ssgsea_kendall$two), invert = as.logical(result_ssgsea_kendall$invert_A))
        result_ssgsea_kendall$P.Value.I=two2one(result_ssgsea_kendall$P.Value, two = as.logical(result_ssgsea_kendall$two), invert = as.logical(result_ssgsea_kendall$invert_I))
        result_ssgsea_kendall_A=result_ssgsea_kendall[,c("geneset_id","P.Value.A")]
        result_ssgsea_kendall_I=result_ssgsea_kendall[,c("geneset_id","P.Value.I")]
        colnames(result_ssgsea_kendall_A)[2]=paste0(vector[i],".ssgsea_kendall")
        colnames(result_ssgsea_kendall_I)[2]=paste0(vector[i],".ssgsea_kendall")
        rm(result_ssgsea_kendall)
      }
      if("ssgsea_spearman" %in% enrichment_method){
        print("ssgsea_spearman will start.")
        result_ssgsea_spearman=spearman_test(geneset=geneset,value=value)
        result_ssgsea_spearman$two=TRUE
        result_ssgsea_spearman$invert_A=ifelse(result_ssgsea_spearman$r<=0,TRUE,FALSE)
        result_ssgsea_spearman$invert_I=ifelse(result_ssgsea_spearman$r>0,TRUE,FALSE)
        result_ssgsea_spearman$P.Value.A=two2one(result_ssgsea_spearman$P.Value, two = as.logical(result_ssgsea_spearman$two), invert = as.logical(result_ssgsea_spearman$invert_A))
        result_ssgsea_spearman$P.Value.I=two2one(result_ssgsea_spearman$P.Value, two = as.logical(result_ssgsea_spearman$two), invert = as.logical(result_ssgsea_spearman$invert_I))
        result_ssgsea_spearman_A=result_ssgsea_spearman[,c("geneset_id","P.Value.A")]
        result_ssgsea_spearman_I=result_ssgsea_spearman[,c("geneset_id","P.Value.I")]
        colnames(result_ssgsea_spearman_A)[2]=paste0(vector[i],".ssgsea_spearman")
        colnames(result_ssgsea_spearman_I)[2]=paste0(vector[i],".ssgsea_spearman")
        rm(result_ssgsea_spearman)
      }
      if("ssgsea_lm" %in% enrichment_method){
        print("ssgsea_lm will start.")
        result_ssgsea_lm=lm_test(geneset=geneset,value=value)
        result_ssgsea_lm$two=TRUE
        result_ssgsea_lm$invert_A=ifelse(result_ssgsea_lm$r<=0,TRUE,FALSE)
        result_ssgsea_lm$invert_I=ifelse(result_ssgsea_lm$r>0,TRUE,FALSE)
        result_ssgsea_lm$P.Value.A=two2one(result_ssgsea_lm$P.Value, two = as.logical(result_ssgsea_lm$two), invert = as.logical(result_ssgsea_lm$invert_A))
        result_ssgsea_lm$P.Value.I=two2one(result_ssgsea_lm$P.Value, two = as.logical(result_ssgsea_lm$two), invert = as.logical(result_ssgsea_lm$invert_I))
        result_ssgsea_lm_A=result_ssgsea_lm[,c("geneset_id","P.Value.A")]
        result_ssgsea_lm_I=result_ssgsea_lm[,c("geneset_id","P.Value.I")]
        colnames(result_ssgsea_lm_A)[2]=paste0(vector[i],".ssgsea_lm")
        colnames(result_ssgsea_lm_I)[2]=paste0(vector[i],".ssgsea_lm")
        rm(result_ssgsea_lm)
      }
      if("ssgsea_biweight" %in% enrichment_method){
        print("ssgsea_biweight will start.")
        result_ssgsea_biweight=biweight_test(geneset=geneset,value=value)
        result_ssgsea_biweight$two=TRUE
        result_ssgsea_biweight$invert_A=ifelse(result_ssgsea_biweight$r<=0,TRUE,FALSE)
        result_ssgsea_biweight$invert_I=ifelse(result_ssgsea_biweight$r>0,TRUE,FALSE)
        result_ssgsea_biweight$P.Value.A=two2one(result_ssgsea_biweight$P.Value, two = as.logical(result_ssgsea_biweight$two), invert = as.logical(result_ssgsea_biweight$invert_A))
        result_ssgsea_biweight$P.Value.I=two2one(result_ssgsea_biweight$P.Value, two = as.logical(result_ssgsea_biweight$two), invert = as.logical(result_ssgsea_biweight$invert_I))
        result_ssgsea_biweight_A=result_ssgsea_biweight[,c("geneset_id","P.Value.A")]
        result_ssgsea_biweight_I=result_ssgsea_biweight[,c("geneset_id","P.Value.I")]
        colnames(result_ssgsea_biweight_A)[2]=paste0(vector[i],".ssgsea_biweight")
        colnames(result_ssgsea_biweight_I)[2]=paste0(vector[i],".ssgsea_biweight")
        rm(result_ssgsea_biweight)
      }
      if("ssgsea_distance" %in% enrichment_method){
        print("ssgsea_distance will start.")
        result_ssgsea_distance=distance_test(geneset=geneset,value=value)
        result_ssgsea_distance$two=TRUE
        result_ssgsea_distance$invert_A=ifelse(result_ssgsea_distance$r<=0,TRUE,FALSE)
        result_ssgsea_distance$invert_I=ifelse(result_ssgsea_distance$r>0,TRUE,FALSE)
        result_ssgsea_distance$P.Value.A=two2one(result_ssgsea_distance$P.Value, two = as.logical(result_ssgsea_distance$two), invert = as.logical(result_ssgsea_distance$invert_A))
        result_ssgsea_distance$P.Value.I=two2one(result_ssgsea_distance$P.Value, two = as.logical(result_ssgsea_distance$two), invert = as.logical(result_ssgsea_distance$invert_I))
        result_ssgsea_distance_A=result_ssgsea_distance[,c("geneset_id","P.Value.A")]
        result_ssgsea_distance_I=result_ssgsea_distance[,c("geneset_id","P.Value.I")]
        colnames(result_ssgsea_distance_A)[2]=paste0(vector[i],".ssgsea_distance")
        colnames(result_ssgsea_distance_I)[2]=paste0(vector[i],".ssgsea_distance")
        rm(result_ssgsea_distance)
      }
      if("ssgsea_percentage" %in% enrichment_method){
        print("ssgsea_percentage will start.")
        result_ssgsea_percentage=percentage_test(geneset=geneset,value=value)
        result_ssgsea_percentage$two=TRUE
        result_ssgsea_percentage$invert_A=ifelse(result_ssgsea_percentage$r<=0,TRUE,FALSE)
        result_ssgsea_percentage$invert_I=ifelse(result_ssgsea_percentage$r>0,TRUE,FALSE)
        result_ssgsea_percentage$P.Value.A=two2one(result_ssgsea_percentage$P.Value, two = as.logical(result_ssgsea_percentage$two), invert = as.logical(result_ssgsea_percentage$invert_A))
        result_ssgsea_percentage$P.Value.I=two2one(result_ssgsea_percentage$P.Value, two = as.logical(result_ssgsea_percentage$two), invert = as.logical(result_ssgsea_percentage$invert_I))
        result_ssgsea_percentage_A=result_ssgsea_percentage[,c("geneset_id","P.Value.A")]
        result_ssgsea_percentage_I=result_ssgsea_percentage[,c("geneset_id","P.Value.I")]
        colnames(result_ssgsea_percentage_A)[2]=paste0(vector[i],".ssgsea_percentage")
        colnames(result_ssgsea_percentage_I)[2]=paste0(vector[i],".ssgsea_percentage")
        rm(result_ssgsea_percentage)
      }
      if("ssgsea_blomqvist" %in% enrichment_method){
        print("ssgsea_blomqvist will start.")
        result_ssgsea_blomqvist=blomqvist_test(geneset=geneset,value=value)
        result_ssgsea_blomqvist$two=TRUE
        result_ssgsea_blomqvist$invert_A=ifelse(result_ssgsea_blomqvist$r<=0,TRUE,FALSE)
        result_ssgsea_blomqvist$invert_I=ifelse(result_ssgsea_blomqvist$r>0,TRUE,FALSE)
        result_ssgsea_blomqvist$P.Value.A=two2one(result_ssgsea_blomqvist$P.Value, two = as.logical(result_ssgsea_blomqvist$two), invert = as.logical(result_ssgsea_blomqvist$invert_A))
        result_ssgsea_blomqvist$P.Value.I=two2one(result_ssgsea_blomqvist$P.Value, two = as.logical(result_ssgsea_blomqvist$two), invert = as.logical(result_ssgsea_blomqvist$invert_I))
        result_ssgsea_blomqvist_A=result_ssgsea_blomqvist[,c("geneset_id","P.Value.A")]
        result_ssgsea_blomqvist_I=result_ssgsea_blomqvist[,c("geneset_id","P.Value.I")]
        colnames(result_ssgsea_blomqvist_A)[2]=paste0(vector[i],".ssgsea_blomqvist")
        colnames(result_ssgsea_blomqvist_I)[2]=paste0(vector[i],".ssgsea_blomqvist")
        rm(result_ssgsea_blomqvist)
      }
      if("ssgsea_hoeffding" %in% enrichment_method){
        print("ssgsea_hoeffding will start.")
        result_ssgsea_hoeffding=hoeffding_test(geneset=geneset,value=value)
        result_ssgsea_hoeffding$two=TRUE
        result_ssgsea_hoeffding$invert_A=ifelse(result_ssgsea_hoeffding$r<=0,TRUE,FALSE)
        result_ssgsea_hoeffding$invert_I=ifelse(result_ssgsea_hoeffding$r>0,TRUE,FALSE)
        result_ssgsea_hoeffding$P.Value.A=two2one(result_ssgsea_hoeffding$P.Value, two = as.logical(result_ssgsea_hoeffding$two), invert = as.logical(result_ssgsea_hoeffding$invert_A))
        result_ssgsea_hoeffding$P.Value.I=two2one(result_ssgsea_hoeffding$P.Value, two = as.logical(result_ssgsea_hoeffding$two), invert = as.logical(result_ssgsea_hoeffding$invert_I))
        result_ssgsea_hoeffding_A=result_ssgsea_hoeffding[,c("geneset_id","P.Value.A")]
        result_ssgsea_hoeffding_I=result_ssgsea_hoeffding[,c("geneset_id","P.Value.I")]
        colnames(result_ssgsea_hoeffding_A)[2]=paste0(vector[i],".ssgsea_hoeffding")
        colnames(result_ssgsea_hoeffding_I)[2]=paste0(vector[i],".ssgsea_hoeffding")
        rm(result_ssgsea_hoeffding)
      }
      if("ssgsea_gamma" %in% enrichment_method){
        print("ssgsea_gamma will start.")
        result_ssgsea_gamma=gamma_test(geneset=geneset,value=value)
        result_ssgsea_gamma$two=TRUE
        result_ssgsea_gamma$invert_A=ifelse(result_ssgsea_gamma$r<=0,TRUE,FALSE)
        result_ssgsea_gamma$invert_I=ifelse(result_ssgsea_gamma$r>0,TRUE,FALSE)
        result_ssgsea_gamma$P.Value.A=two2one(result_ssgsea_gamma$P.Value, two = as.logical(result_ssgsea_gamma$two), invert = as.logical(result_ssgsea_gamma$invert_A))
        result_ssgsea_gamma$P.Value.I=two2one(result_ssgsea_gamma$P.Value, two = as.logical(result_ssgsea_gamma$two), invert = as.logical(result_ssgsea_gamma$invert_I))
        result_ssgsea_gamma_A=result_ssgsea_gamma[,c("geneset_id","P.Value.A")]
        result_ssgsea_gamma_I=result_ssgsea_gamma[,c("geneset_id","P.Value.I")]
        colnames(result_ssgsea_gamma_A)[2]=paste0(vector[i],".ssgsea_gamma")
        colnames(result_ssgsea_gamma_I)[2]=paste0(vector[i],".ssgsea_gamma")
        rm(result_ssgsea_gamma)
      }
    }else {
      print("Not using ssgsea.")
    }
    ##########################################zscore
    methods_to_check=c("zscore_pearson","zscore_kendall","zscore_spearman","zscore_lm","zscore_biweight","zscore_distance","zscore_percentage","zscore_blomqvist","zscore_hoeffding","zscore_gamma")
    if (any(methods_to_check %in% enrichment_method)){
      geneset=zscoreParam(exprData=as.matrix(exp),geneSets=geneSets,minSize=min.sz,maxSize=max.sz)
      geneset=gsva(geneset)
      geneset=as.data.frame(geneset)
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("zscore_pearson" %in% enrichment_method){
        print("zscore_pearson will start.")
        result_zscore_pearson=pearson_test(geneset=geneset,value=value)
        result_zscore_pearson$two=TRUE
        result_zscore_pearson$invert_A=ifelse(result_zscore_pearson$r<=0,TRUE,FALSE)
        result_zscore_pearson$invert_I=ifelse(result_zscore_pearson$r>0,TRUE,FALSE)
        result_zscore_pearson$P.Value.A=two2one(result_zscore_pearson$P.Value, two = as.logical(result_zscore_pearson$two), invert = as.logical(result_zscore_pearson$invert_A))
        result_zscore_pearson$P.Value.I=two2one(result_zscore_pearson$P.Value, two = as.logical(result_zscore_pearson$two), invert = as.logical(result_zscore_pearson$invert_I))
        result_zscore_pearson_A=result_zscore_pearson[,c("geneset_id","P.Value.A")]
        result_zscore_pearson_I=result_zscore_pearson[,c("geneset_id","P.Value.I")]
        colnames(result_zscore_pearson_A)[2]=paste0(vector[i],".zscore_pearson")
        colnames(result_zscore_pearson_I)[2]=paste0(vector[i],".zscore_pearson")
        rm(result_zscore_pearson)
      }
      if("zscore_kendall" %in% enrichment_method){
        print("zscore_kendall will start.")
        result_zscore_kendall=kendall_test(geneset=geneset,value=value)
        result_zscore_kendall$two=TRUE
        result_zscore_kendall$invert_A=ifelse(result_zscore_kendall$r<=0,TRUE,FALSE)
        result_zscore_kendall$invert_I=ifelse(result_zscore_kendall$r>0,TRUE,FALSE)
        result_zscore_kendall$P.Value.A=two2one(result_zscore_kendall$P.Value, two = as.logical(result_zscore_kendall$two), invert = as.logical(result_zscore_kendall$invert_A))
        result_zscore_kendall$P.Value.I=two2one(result_zscore_kendall$P.Value, two = as.logical(result_zscore_kendall$two), invert = as.logical(result_zscore_kendall$invert_I))
        result_zscore_kendall_A=result_zscore_kendall[,c("geneset_id","P.Value.A")]
        result_zscore_kendall_I=result_zscore_kendall[,c("geneset_id","P.Value.I")]
        colnames(result_zscore_kendall_A)[2]=paste0(vector[i],".zscore_kendall")
        colnames(result_zscore_kendall_I)[2]=paste0(vector[i],".zscore_kendall")
        rm(result_zscore_kendall)
      }
      if("zscore_spearman" %in% enrichment_method){
        print("zscore_spearman will start.")
        result_zscore_spearman=spearman_test(geneset=geneset,value=value)
        result_zscore_spearman$two=TRUE
        result_zscore_spearman$invert_A=ifelse(result_zscore_spearman$r<=0,TRUE,FALSE)
        result_zscore_spearman$invert_I=ifelse(result_zscore_spearman$r>0,TRUE,FALSE)
        result_zscore_spearman$P.Value.A=two2one(result_zscore_spearman$P.Value, two = as.logical(result_zscore_spearman$two), invert = as.logical(result_zscore_spearman$invert_A))
        result_zscore_spearman$P.Value.I=two2one(result_zscore_spearman$P.Value, two = as.logical(result_zscore_spearman$two), invert = as.logical(result_zscore_spearman$invert_I))
        result_zscore_spearman_A=result_zscore_spearman[,c("geneset_id","P.Value.A")]
        result_zscore_spearman_I=result_zscore_spearman[,c("geneset_id","P.Value.I")]
        colnames(result_zscore_spearman_A)[2]=paste0(vector[i],".zscore_spearman")
        colnames(result_zscore_spearman_I)[2]=paste0(vector[i],".zscore_spearman")
        rm(result_zscore_spearman)
      }
      if("zscore_lm" %in% enrichment_method){
        print("zscore_lm will start.")
        result_zscore_lm=lm_test(geneset=geneset,value=value)
        result_zscore_lm$two=TRUE
        result_zscore_lm$invert_A=ifelse(result_zscore_lm$r<=0,TRUE,FALSE)
        result_zscore_lm$invert_I=ifelse(result_zscore_lm$r>0,TRUE,FALSE)
        result_zscore_lm$P.Value.A=two2one(result_zscore_lm$P.Value, two = as.logical(result_zscore_lm$two), invert = as.logical(result_zscore_lm$invert_A))
        result_zscore_lm$P.Value.I=two2one(result_zscore_lm$P.Value, two = as.logical(result_zscore_lm$two), invert = as.logical(result_zscore_lm$invert_I))
        result_zscore_lm_A=result_zscore_lm[,c("geneset_id","P.Value.A")]
        result_zscore_lm_I=result_zscore_lm[,c("geneset_id","P.Value.I")]
        colnames(result_zscore_lm_A)[2]=paste0(vector[i],".zscore_lm")
        colnames(result_zscore_lm_I)[2]=paste0(vector[i],".zscore_lm")
        rm(result_zscore_lm)
      }
      if("zscore_biweight" %in% enrichment_method){
        print("zscore_biweight will start.")
        result_zscore_biweight=biweight_test(geneset=geneset,value=value)
        result_zscore_biweight$two=TRUE
        result_zscore_biweight$invert_A=ifelse(result_zscore_biweight$r<=0,TRUE,FALSE)
        result_zscore_biweight$invert_I=ifelse(result_zscore_biweight$r>0,TRUE,FALSE)
        result_zscore_biweight$P.Value.A=two2one(result_zscore_biweight$P.Value, two = as.logical(result_zscore_biweight$two), invert = as.logical(result_zscore_biweight$invert_A))
        result_zscore_biweight$P.Value.I=two2one(result_zscore_biweight$P.Value, two = as.logical(result_zscore_biweight$two), invert = as.logical(result_zscore_biweight$invert_I))
        result_zscore_biweight_A=result_zscore_biweight[,c("geneset_id","P.Value.A")]
        result_zscore_biweight_I=result_zscore_biweight[,c("geneset_id","P.Value.I")]
        colnames(result_zscore_biweight_A)[2]=paste0(vector[i],".zscore_biweight")
        colnames(result_zscore_biweight_I)[2]=paste0(vector[i],".zscore_biweight")
        rm(result_zscore_biweight)
      }
      if("zscore_distance" %in% enrichment_method){
        print("zscore_distance will start.")
        result_zscore_distance=distance_test(geneset=geneset,value=value)
        result_zscore_distance$two=TRUE
        result_zscore_distance$invert_A=ifelse(result_zscore_distance$r<=0,TRUE,FALSE)
        result_zscore_distance$invert_I=ifelse(result_zscore_distance$r>0,TRUE,FALSE)
        result_zscore_distance$P.Value.A=two2one(result_zscore_distance$P.Value, two = as.logical(result_zscore_distance$two), invert = as.logical(result_zscore_distance$invert_A))
        result_zscore_distance$P.Value.I=two2one(result_zscore_distance$P.Value, two = as.logical(result_zscore_distance$two), invert = as.logical(result_zscore_distance$invert_I))
        result_zscore_distance_A=result_zscore_distance[,c("geneset_id","P.Value.A")]
        result_zscore_distance_I=result_zscore_distance[,c("geneset_id","P.Value.I")]
        colnames(result_zscore_distance_A)[2]=paste0(vector[i],".zscore_distance")
        colnames(result_zscore_distance_I)[2]=paste0(vector[i],".zscore_distance")
        rm(result_zscore_distance)
      }
      if("zscore_percentage" %in% enrichment_method){
        print("zscore_percentage will start.")
        result_zscore_percentage=percentage_test(geneset=geneset,value=value)
        result_zscore_percentage$two=TRUE
        result_zscore_percentage$invert_A=ifelse(result_zscore_percentage$r<=0,TRUE,FALSE)
        result_zscore_percentage$invert_I=ifelse(result_zscore_percentage$r>0,TRUE,FALSE)
        result_zscore_percentage$P.Value.A=two2one(result_zscore_percentage$P.Value, two = as.logical(result_zscore_percentage$two), invert = as.logical(result_zscore_percentage$invert_A))
        result_zscore_percentage$P.Value.I=two2one(result_zscore_percentage$P.Value, two = as.logical(result_zscore_percentage$two), invert = as.logical(result_zscore_percentage$invert_I))
        result_zscore_percentage_A=result_zscore_percentage[,c("geneset_id","P.Value.A")]
        result_zscore_percentage_I=result_zscore_percentage[,c("geneset_id","P.Value.I")]
        colnames(result_zscore_percentage_A)[2]=paste0(vector[i],".zscore_percentage")
        colnames(result_zscore_percentage_I)[2]=paste0(vector[i],".zscore_percentage")
        rm(result_zscore_percentage)
      }
      if("zscore_blomqvist" %in% enrichment_method){
        print("zscore_blomqvist will start.")
        result_zscore_blomqvist=blomqvist_test(geneset=geneset,value=value)
        result_zscore_blomqvist$two=TRUE
        result_zscore_blomqvist$invert_A=ifelse(result_zscore_blomqvist$r<=0,TRUE,FALSE)
        result_zscore_blomqvist$invert_I=ifelse(result_zscore_blomqvist$r>0,TRUE,FALSE)
        result_zscore_blomqvist$P.Value.A=two2one(result_zscore_blomqvist$P.Value, two = as.logical(result_zscore_blomqvist$two), invert = as.logical(result_zscore_blomqvist$invert_A))
        result_zscore_blomqvist$P.Value.I=two2one(result_zscore_blomqvist$P.Value, two = as.logical(result_zscore_blomqvist$two), invert = as.logical(result_zscore_blomqvist$invert_I))
        result_zscore_blomqvist_A=result_zscore_blomqvist[,c("geneset_id","P.Value.A")]
        result_zscore_blomqvist_I=result_zscore_blomqvist[,c("geneset_id","P.Value.I")]
        colnames(result_zscore_blomqvist_A)[2]=paste0(vector[i],".zscore_blomqvist")
        colnames(result_zscore_blomqvist_I)[2]=paste0(vector[i],".zscore_blomqvist")
        rm(result_zscore_blomqvist)
      }
      if("zscore_hoeffding" %in% enrichment_method){
        print("zscore_hoeffding will start.")
        result_zscore_hoeffding=hoeffding_test(geneset=geneset,value=value)
        result_zscore_hoeffding$two=TRUE
        result_zscore_hoeffding$invert_A=ifelse(result_zscore_hoeffding$r<=0,TRUE,FALSE)
        result_zscore_hoeffding$invert_I=ifelse(result_zscore_hoeffding$r>0,TRUE,FALSE)
        result_zscore_hoeffding$P.Value.A=two2one(result_zscore_hoeffding$P.Value, two = as.logical(result_zscore_hoeffding$two), invert = as.logical(result_zscore_hoeffding$invert_A))
        result_zscore_hoeffding$P.Value.I=two2one(result_zscore_hoeffding$P.Value, two = as.logical(result_zscore_hoeffding$two), invert = as.logical(result_zscore_hoeffding$invert_I))
        result_zscore_hoeffding_A=result_zscore_hoeffding[,c("geneset_id","P.Value.A")]
        result_zscore_hoeffding_I=result_zscore_hoeffding[,c("geneset_id","P.Value.I")]
        colnames(result_zscore_hoeffding_A)[2]=paste0(vector[i],".zscore_hoeffding")
        colnames(result_zscore_hoeffding_I)[2]=paste0(vector[i],".zscore_hoeffding")
        rm(result_zscore_hoeffding)
      }
      if("zscore_gamma" %in% enrichment_method){
        print("zscore_gamma will start.")
        result_zscore_gamma=gamma_test(geneset=geneset,value=value)
        result_zscore_gamma$two=TRUE
        result_zscore_gamma$invert_A=ifelse(result_zscore_gamma$r<=0,TRUE,FALSE)
        result_zscore_gamma$invert_I=ifelse(result_zscore_gamma$r>0,TRUE,FALSE)
        result_zscore_gamma$P.Value.A=two2one(result_zscore_gamma$P.Value, two = as.logical(result_zscore_gamma$two), invert = as.logical(result_zscore_gamma$invert_A))
        result_zscore_gamma$P.Value.I=two2one(result_zscore_gamma$P.Value, two = as.logical(result_zscore_gamma$two), invert = as.logical(result_zscore_gamma$invert_I))
        result_zscore_gamma_A=result_zscore_gamma[,c("geneset_id","P.Value.A")]
        result_zscore_gamma_I=result_zscore_gamma[,c("geneset_id","P.Value.I")]
        colnames(result_zscore_gamma_A)[2]=paste0(vector[i],".zscore_gamma")
        colnames(result_zscore_gamma_I)[2]=paste0(vector[i],".zscore_gamma")
        rm(result_zscore_gamma)
      }
    }else {
      print("Not using zscore.")
    }
    ##########################################plage
    methods_to_check=c("plage_pearson","plage_kendall","plage_spearman","plage_lm","plage_biweight","plage_distance","plage_percentage","plage_blomqvist","plage_hoeffding","plage_gamma")
    if (any(methods_to_check %in% enrichment_method)){
      geneset=plageParam(exprData=as.matrix(exp),geneSets=geneSets,minSize=min.sz,maxSize=max.sz)
      geneset=gsva(geneset)
      geneset=as.data.frame(geneset)
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("plage_pearson" %in% enrichment_method){
        print("plage_pearson will start.")
        result_plage_pearson=pearson_test(geneset=geneset,value=value)
        result_plage_pearson$two=TRUE
        result_plage_pearson$invert_A=ifelse(result_plage_pearson$r<=0,TRUE,FALSE)
        result_plage_pearson$invert_I=ifelse(result_plage_pearson$r>0,TRUE,FALSE)
        result_plage_pearson$P.Value.A=two2one(result_plage_pearson$P.Value, two = as.logical(result_plage_pearson$two), invert = as.logical(result_plage_pearson$invert_A))
        result_plage_pearson$P.Value.I=two2one(result_plage_pearson$P.Value, two = as.logical(result_plage_pearson$two), invert = as.logical(result_plage_pearson$invert_I))
        result_plage_pearson_A=result_plage_pearson[,c("geneset_id","P.Value.A")]
        result_plage_pearson_I=result_plage_pearson[,c("geneset_id","P.Value.I")]
        colnames(result_plage_pearson_A)[2]=paste0(vector[i],".plage_pearson")
        colnames(result_plage_pearson_I)[2]=paste0(vector[i],".plage_pearson")
        rm(result_plage_pearson)
      }
      if("plage_kendall" %in% enrichment_method){
        print("plage_kendall will start.")
        result_plage_kendall=kendall_test(geneset=geneset,value=value)
        result_plage_kendall$two=TRUE
        result_plage_kendall$invert_A=ifelse(result_plage_kendall$r<=0,TRUE,FALSE)
        result_plage_kendall$invert_I=ifelse(result_plage_kendall$r>0,TRUE,FALSE)
        result_plage_kendall$P.Value.A=two2one(result_plage_kendall$P.Value, two = as.logical(result_plage_kendall$two), invert = as.logical(result_plage_kendall$invert_A))
        result_plage_kendall$P.Value.I=two2one(result_plage_kendall$P.Value, two = as.logical(result_plage_kendall$two), invert = as.logical(result_plage_kendall$invert_I))
        result_plage_kendall_A=result_plage_kendall[,c("geneset_id","P.Value.A")]
        result_plage_kendall_I=result_plage_kendall[,c("geneset_id","P.Value.I")]
        colnames(result_plage_kendall_A)[2]=paste0(vector[i],".plage_kendall")
        colnames(result_plage_kendall_I)[2]=paste0(vector[i],".plage_kendall")
        rm(result_plage_kendall)
      }
      if("plage_spearman" %in% enrichment_method){
        print("plage_spearman will start.")
        result_plage_spearman=spearman_test(geneset=geneset,value=value)
        result_plage_spearman$two=TRUE
        result_plage_spearman$invert_A=ifelse(result_plage_spearman$r<=0,TRUE,FALSE)
        result_plage_spearman$invert_I=ifelse(result_plage_spearman$r>0,TRUE,FALSE)
        result_plage_spearman$P.Value.A=two2one(result_plage_spearman$P.Value, two = as.logical(result_plage_spearman$two), invert = as.logical(result_plage_spearman$invert_A))
        result_plage_spearman$P.Value.I=two2one(result_plage_spearman$P.Value, two = as.logical(result_plage_spearman$two), invert = as.logical(result_plage_spearman$invert_I))
        result_plage_spearman_A=result_plage_spearman[,c("geneset_id","P.Value.A")]
        result_plage_spearman_I=result_plage_spearman[,c("geneset_id","P.Value.I")]
        colnames(result_plage_spearman_A)[2]=paste0(vector[i],".plage_spearman")
        colnames(result_plage_spearman_I)[2]=paste0(vector[i],".plage_spearman")
        rm(result_plage_spearman)
      }
      if("plage_lm" %in% enrichment_method){
        print("plage_lm will start.")
        result_plage_lm=lm_test(geneset=geneset,value=value)
        result_plage_lm$two=TRUE
        result_plage_lm$invert_A=ifelse(result_plage_lm$r<=0,TRUE,FALSE)
        result_plage_lm$invert_I=ifelse(result_plage_lm$r>0,TRUE,FALSE)
        result_plage_lm$P.Value.A=two2one(result_plage_lm$P.Value, two = as.logical(result_plage_lm$two), invert = as.logical(result_plage_lm$invert_A))
        result_plage_lm$P.Value.I=two2one(result_plage_lm$P.Value, two = as.logical(result_plage_lm$two), invert = as.logical(result_plage_lm$invert_I))
        result_plage_lm_A=result_plage_lm[,c("geneset_id","P.Value.A")]
        result_plage_lm_I=result_plage_lm[,c("geneset_id","P.Value.I")]
        colnames(result_plage_lm_A)[2]=paste0(vector[i],".plage_lm")
        colnames(result_plage_lm_I)[2]=paste0(vector[i],".plage_lm")
        rm(result_plage_lm)
      }
      if("plage_biweight" %in% enrichment_method){
        print("plage_biweight will start.")
        result_plage_biweight=biweight_test(geneset=geneset,value=value)
        result_plage_biweight$two=TRUE
        result_plage_biweight$invert_A=ifelse(result_plage_biweight$r<=0,TRUE,FALSE)
        result_plage_biweight$invert_I=ifelse(result_plage_biweight$r>0,TRUE,FALSE)
        result_plage_biweight$P.Value.A=two2one(result_plage_biweight$P.Value, two = as.logical(result_plage_biweight$two), invert = as.logical(result_plage_biweight$invert_A))
        result_plage_biweight$P.Value.I=two2one(result_plage_biweight$P.Value, two = as.logical(result_plage_biweight$two), invert = as.logical(result_plage_biweight$invert_I))
        result_plage_biweight_A=result_plage_biweight[,c("geneset_id","P.Value.A")]
        result_plage_biweight_I=result_plage_biweight[,c("geneset_id","P.Value.I")]
        colnames(result_plage_biweight_A)[2]=paste0(vector[i],".plage_biweight")
        colnames(result_plage_biweight_I)[2]=paste0(vector[i],".plage_biweight")
        rm(result_plage_biweight)
      }
      if("plage_distance" %in% enrichment_method){
        print("plage_distance will start.")
        result_plage_distance=distance_test(geneset=geneset,value=value)
        result_plage_distance$two=TRUE
        result_plage_distance$invert_A=ifelse(result_plage_distance$r<=0,TRUE,FALSE)
        result_plage_distance$invert_I=ifelse(result_plage_distance$r>0,TRUE,FALSE)
        result_plage_distance$P.Value.A=two2one(result_plage_distance$P.Value, two = as.logical(result_plage_distance$two), invert = as.logical(result_plage_distance$invert_A))
        result_plage_distance$P.Value.I=two2one(result_plage_distance$P.Value, two = as.logical(result_plage_distance$two), invert = as.logical(result_plage_distance$invert_I))
        result_plage_distance_A=result_plage_distance[,c("geneset_id","P.Value.A")]
        result_plage_distance_I=result_plage_distance[,c("geneset_id","P.Value.I")]
        colnames(result_plage_distance_A)[2]=paste0(vector[i],".plage_distance")
        colnames(result_plage_distance_I)[2]=paste0(vector[i],".plage_distance")
        rm(result_plage_distance)
      }
      if("plage_percentage" %in% enrichment_method){
        print("plage_percentage will start.")
        result_plage_percentage=percentage_test(geneset=geneset,value=value)
        result_plage_percentage$two=TRUE
        result_plage_percentage$invert_A=ifelse(result_plage_percentage$r<=0,TRUE,FALSE)
        result_plage_percentage$invert_I=ifelse(result_plage_percentage$r>0,TRUE,FALSE)
        result_plage_percentage$P.Value.A=two2one(result_plage_percentage$P.Value, two = as.logical(result_plage_percentage$two), invert = as.logical(result_plage_percentage$invert_A))
        result_plage_percentage$P.Value.I=two2one(result_plage_percentage$P.Value, two = as.logical(result_plage_percentage$two), invert = as.logical(result_plage_percentage$invert_I))
        result_plage_percentage_A=result_plage_percentage[,c("geneset_id","P.Value.A")]
        result_plage_percentage_I=result_plage_percentage[,c("geneset_id","P.Value.I")]
        colnames(result_plage_percentage_A)[2]=paste0(vector[i],".plage_percentage")
        colnames(result_plage_percentage_I)[2]=paste0(vector[i],".plage_percentage")
        rm(result_plage_percentage)
      }
      if("plage_blomqvist" %in% enrichment_method){
        print("plage_blomqvist will start.")
        result_plage_blomqvist=blomqvist_test(geneset=geneset,value=value)
        result_plage_blomqvist$two=TRUE
        result_plage_blomqvist$invert_A=ifelse(result_plage_blomqvist$r<=0,TRUE,FALSE)
        result_plage_blomqvist$invert_I=ifelse(result_plage_blomqvist$r>0,TRUE,FALSE)
        result_plage_blomqvist$P.Value.A=two2one(result_plage_blomqvist$P.Value, two = as.logical(result_plage_blomqvist$two), invert = as.logical(result_plage_blomqvist$invert_A))
        result_plage_blomqvist$P.Value.I=two2one(result_plage_blomqvist$P.Value, two = as.logical(result_plage_blomqvist$two), invert = as.logical(result_plage_blomqvist$invert_I))
        result_plage_blomqvist_A=result_plage_blomqvist[,c("geneset_id","P.Value.A")]
        result_plage_blomqvist_I=result_plage_blomqvist[,c("geneset_id","P.Value.I")]
        colnames(result_plage_blomqvist_A)[2]=paste0(vector[i],".plage_blomqvist")
        colnames(result_plage_blomqvist_I)[2]=paste0(vector[i],".plage_blomqvist")
        rm(result_plage_blomqvist)
      }
      if("plage_hoeffding" %in% enrichment_method){
        print("plage_hoeffding will start.")
        result_plage_hoeffding=hoeffding_test(geneset=geneset,value=value)
        result_plage_hoeffding$two=TRUE
        result_plage_hoeffding$invert_A=ifelse(result_plage_hoeffding$r<=0,TRUE,FALSE)
        result_plage_hoeffding$invert_I=ifelse(result_plage_hoeffding$r>0,TRUE,FALSE)
        result_plage_hoeffding$P.Value.A=two2one(result_plage_hoeffding$P.Value, two = as.logical(result_plage_hoeffding$two), invert = as.logical(result_plage_hoeffding$invert_A))
        result_plage_hoeffding$P.Value.I=two2one(result_plage_hoeffding$P.Value, two = as.logical(result_plage_hoeffding$two), invert = as.logical(result_plage_hoeffding$invert_I))
        result_plage_hoeffding_A=result_plage_hoeffding[,c("geneset_id","P.Value.A")]
        result_plage_hoeffding_I=result_plage_hoeffding[,c("geneset_id","P.Value.I")]
        colnames(result_plage_hoeffding_A)[2]=paste0(vector[i],".plage_hoeffding")
        colnames(result_plage_hoeffding_I)[2]=paste0(vector[i],".plage_hoeffding")
        rm(result_plage_hoeffding)
      }
      if("plage_gamma" %in% enrichment_method){
        print("plage_gamma will start.")
        result_plage_gamma=gamma_test(geneset=geneset,value=value)
        result_plage_gamma$two=TRUE
        result_plage_gamma$invert_A=ifelse(result_plage_gamma$r<=0,TRUE,FALSE)
        result_plage_gamma$invert_I=ifelse(result_plage_gamma$r>0,TRUE,FALSE)
        result_plage_gamma$P.Value.A=two2one(result_plage_gamma$P.Value, two = as.logical(result_plage_gamma$two), invert = as.logical(result_plage_gamma$invert_A))
        result_plage_gamma$P.Value.I=two2one(result_plage_gamma$P.Value, two = as.logical(result_plage_gamma$two), invert = as.logical(result_plage_gamma$invert_I))
        result_plage_gamma_A=result_plage_gamma[,c("geneset_id","P.Value.A")]
        result_plage_gamma_I=result_plage_gamma[,c("geneset_id","P.Value.I")]
        colnames(result_plage_gamma_A)[2]=paste0(vector[i],".plage_gamma")
        colnames(result_plage_gamma_I)[2]=paste0(vector[i],".plage_gamma")
        rm(result_plage_gamma)
      }
    }else {
      print("Not using plage.")
    }
    ##########################################pca
    methods_to_check=c("pca_pearson","pca_kendall","pca_spearman","pca_lm","pca_biweight","pca_distance","pca_percentage","pca_blomqvist","pca_hoeffding","pca_gamma")
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
      if("pca_pearson" %in% enrichment_method){
        print("pca_pearson will start.")
        result_pca_pearson=pearson_test(geneset=geneset,value=value)
        result_pca_pearson$two=TRUE
        result_pca_pearson$invert_A=ifelse(result_pca_pearson$r<=0,TRUE,FALSE)
        result_pca_pearson$invert_I=ifelse(result_pca_pearson$r>0,TRUE,FALSE)
        result_pca_pearson$P.Value.A=two2one(result_pca_pearson$P.Value, two = as.logical(result_pca_pearson$two), invert = as.logical(result_pca_pearson$invert_A))
        result_pca_pearson$P.Value.I=two2one(result_pca_pearson$P.Value, two = as.logical(result_pca_pearson$two), invert = as.logical(result_pca_pearson$invert_I))
        result_pca_pearson_A=result_pca_pearson[,c("geneset_id","P.Value.A")]
        result_pca_pearson_I=result_pca_pearson[,c("geneset_id","P.Value.I")]
        colnames(result_pca_pearson_A)[2]=paste0(vector[i],".pca_pearson")
        colnames(result_pca_pearson_I)[2]=paste0(vector[i],".pca_pearson")
        rm(result_pca_pearson)
      }
      if("pca_kendall" %in% enrichment_method){
        print("pca_kendall will start.")
        result_pca_kendall=kendall_test(geneset=geneset,value=value)
        result_pca_kendall$two=TRUE
        result_pca_kendall$invert_A=ifelse(result_pca_kendall$r<=0,TRUE,FALSE)
        result_pca_kendall$invert_I=ifelse(result_pca_kendall$r>0,TRUE,FALSE)
        result_pca_kendall$P.Value.A=two2one(result_pca_kendall$P.Value, two = as.logical(result_pca_kendall$two), invert = as.logical(result_pca_kendall$invert_A))
        result_pca_kendall$P.Value.I=two2one(result_pca_kendall$P.Value, two = as.logical(result_pca_kendall$two), invert = as.logical(result_pca_kendall$invert_I))
        result_pca_kendall_A=result_pca_kendall[,c("geneset_id","P.Value.A")]
        result_pca_kendall_I=result_pca_kendall[,c("geneset_id","P.Value.I")]
        colnames(result_pca_kendall_A)[2]=paste0(vector[i],".pca_kendall")
        colnames(result_pca_kendall_I)[2]=paste0(vector[i],".pca_kendall")
        rm(result_pca_kendall)
      }
      if("pca_spearman" %in% enrichment_method){
        print("pca_spearman will start.")
        result_pca_spearman=spearman_test(geneset=geneset,value=value)
        result_pca_spearman$two=TRUE
        result_pca_spearman$invert_A=ifelse(result_pca_spearman$r<=0,TRUE,FALSE)
        result_pca_spearman$invert_I=ifelse(result_pca_spearman$r>0,TRUE,FALSE)
        result_pca_spearman$P.Value.A=two2one(result_pca_spearman$P.Value, two = as.logical(result_pca_spearman$two), invert = as.logical(result_pca_spearman$invert_A))
        result_pca_spearman$P.Value.I=two2one(result_pca_spearman$P.Value, two = as.logical(result_pca_spearman$two), invert = as.logical(result_pca_spearman$invert_I))
        result_pca_spearman_A=result_pca_spearman[,c("geneset_id","P.Value.A")]
        result_pca_spearman_I=result_pca_spearman[,c("geneset_id","P.Value.I")]
        colnames(result_pca_spearman_A)[2]=paste0(vector[i],".pca_spearman")
        colnames(result_pca_spearman_I)[2]=paste0(vector[i],".pca_spearman")
        rm(result_pca_spearman)
      }
      if("pca_lm" %in% enrichment_method){
        print("pca_lm will start.")
        result_pca_lm=lm_test(geneset=geneset,value=value)
        result_pca_lm$two=TRUE
        result_pca_lm$invert_A=ifelse(result_pca_lm$r<=0,TRUE,FALSE)
        result_pca_lm$invert_I=ifelse(result_pca_lm$r>0,TRUE,FALSE)
        result_pca_lm$P.Value.A=two2one(result_pca_lm$P.Value, two = as.logical(result_pca_lm$two), invert = as.logical(result_pca_lm$invert_A))
        result_pca_lm$P.Value.I=two2one(result_pca_lm$P.Value, two = as.logical(result_pca_lm$two), invert = as.logical(result_pca_lm$invert_I))
        result_pca_lm_A=result_pca_lm[,c("geneset_id","P.Value.A")]
        result_pca_lm_I=result_pca_lm[,c("geneset_id","P.Value.I")]
        colnames(result_pca_lm_A)[2]=paste0(vector[i],".pca_lm")
        colnames(result_pca_lm_I)[2]=paste0(vector[i],".pca_lm")
        rm(result_pca_lm)
      }
      if("pca_biweight" %in% enrichment_method){
        print("pca_biweight will start.")
        result_pca_biweight=biweight_test(geneset=geneset,value=value)
        result_pca_biweight$two=TRUE
        result_pca_biweight$invert_A=ifelse(result_pca_biweight$r<=0,TRUE,FALSE)
        result_pca_biweight$invert_I=ifelse(result_pca_biweight$r>0,TRUE,FALSE)
        result_pca_biweight$P.Value.A=two2one(result_pca_biweight$P.Value, two = as.logical(result_pca_biweight$two), invert = as.logical(result_pca_biweight$invert_A))
        result_pca_biweight$P.Value.I=two2one(result_pca_biweight$P.Value, two = as.logical(result_pca_biweight$two), invert = as.logical(result_pca_biweight$invert_I))
        result_pca_biweight_A=result_pca_biweight[,c("geneset_id","P.Value.A")]
        result_pca_biweight_I=result_pca_biweight[,c("geneset_id","P.Value.I")]
        colnames(result_pca_biweight_A)[2]=paste0(vector[i],".pca_biweight")
        colnames(result_pca_biweight_I)[2]=paste0(vector[i],".pca_biweight")
        rm(result_pca_biweight)
      }
      if("pca_distance" %in% enrichment_method){
        print("pca_distance will start.")
        result_pca_distance=distance_test(geneset=geneset,value=value)
        result_pca_distance$two=TRUE
        result_pca_distance$invert_A=ifelse(result_pca_distance$r<=0,TRUE,FALSE)
        result_pca_distance$invert_I=ifelse(result_pca_distance$r>0,TRUE,FALSE)
        result_pca_distance$P.Value.A=two2one(result_pca_distance$P.Value, two = as.logical(result_pca_distance$two), invert = as.logical(result_pca_distance$invert_A))
        result_pca_distance$P.Value.I=two2one(result_pca_distance$P.Value, two = as.logical(result_pca_distance$two), invert = as.logical(result_pca_distance$invert_I))
        result_pca_distance_A=result_pca_distance[,c("geneset_id","P.Value.A")]
        result_pca_distance_I=result_pca_distance[,c("geneset_id","P.Value.I")]
        colnames(result_pca_distance_A)[2]=paste0(vector[i],".pca_distance")
        colnames(result_pca_distance_I)[2]=paste0(vector[i],".pca_distance")
        rm(result_pca_distance)
      }
      if("pca_percentage" %in% enrichment_method){
        print("pca_percentage will start.")
        result_pca_percentage=percentage_test(geneset=geneset,value=value)
        result_pca_percentage$two=TRUE
        result_pca_percentage$invert_A=ifelse(result_pca_percentage$r<=0,TRUE,FALSE)
        result_pca_percentage$invert_I=ifelse(result_pca_percentage$r>0,TRUE,FALSE)
        result_pca_percentage$P.Value.A=two2one(result_pca_percentage$P.Value, two = as.logical(result_pca_percentage$two), invert = as.logical(result_pca_percentage$invert_A))
        result_pca_percentage$P.Value.I=two2one(result_pca_percentage$P.Value, two = as.logical(result_pca_percentage$two), invert = as.logical(result_pca_percentage$invert_I))
        result_pca_percentage_A=result_pca_percentage[,c("geneset_id","P.Value.A")]
        result_pca_percentage_I=result_pca_percentage[,c("geneset_id","P.Value.I")]
        colnames(result_pca_percentage_A)[2]=paste0(vector[i],".pca_percentage")
        colnames(result_pca_percentage_I)[2]=paste0(vector[i],".pca_percentage")
        rm(result_pca_percentage)
      }
      if("pca_blomqvist" %in% enrichment_method){
        print("pca_blomqvist will start.")
        result_pca_blomqvist=blomqvist_test(geneset=geneset,value=value)
        result_pca_blomqvist$two=TRUE
        result_pca_blomqvist$invert_A=ifelse(result_pca_blomqvist$r<=0,TRUE,FALSE)
        result_pca_blomqvist$invert_I=ifelse(result_pca_blomqvist$r>0,TRUE,FALSE)
        result_pca_blomqvist$P.Value.A=two2one(result_pca_blomqvist$P.Value, two = as.logical(result_pca_blomqvist$two), invert = as.logical(result_pca_blomqvist$invert_A))
        result_pca_blomqvist$P.Value.I=two2one(result_pca_blomqvist$P.Value, two = as.logical(result_pca_blomqvist$two), invert = as.logical(result_pca_blomqvist$invert_I))
        result_pca_blomqvist_A=result_pca_blomqvist[,c("geneset_id","P.Value.A")]
        result_pca_blomqvist_I=result_pca_blomqvist[,c("geneset_id","P.Value.I")]
        colnames(result_pca_blomqvist_A)[2]=paste0(vector[i],".pca_blomqvist")
        colnames(result_pca_blomqvist_I)[2]=paste0(vector[i],".pca_blomqvist")
        rm(result_pca_blomqvist)
      }
      if("pca_hoeffding" %in% enrichment_method){
        print("pca_hoeffding will start.")
        result_pca_hoeffding=hoeffding_test(geneset=geneset,value=value)
        result_pca_hoeffding$two=TRUE
        result_pca_hoeffding$invert_A=ifelse(result_pca_hoeffding$r<=0,TRUE,FALSE)
        result_pca_hoeffding$invert_I=ifelse(result_pca_hoeffding$r>0,TRUE,FALSE)
        result_pca_hoeffding$P.Value.A=two2one(result_pca_hoeffding$P.Value, two = as.logical(result_pca_hoeffding$two), invert = as.logical(result_pca_hoeffding$invert_A))
        result_pca_hoeffding$P.Value.I=two2one(result_pca_hoeffding$P.Value, two = as.logical(result_pca_hoeffding$two), invert = as.logical(result_pca_hoeffding$invert_I))
        result_pca_hoeffding_A=result_pca_hoeffding[,c("geneset_id","P.Value.A")]
        result_pca_hoeffding_I=result_pca_hoeffding[,c("geneset_id","P.Value.I")]
        colnames(result_pca_hoeffding_A)[2]=paste0(vector[i],".pca_hoeffding")
        colnames(result_pca_hoeffding_I)[2]=paste0(vector[i],".pca_hoeffding")
        rm(result_pca_hoeffding)
      }
      if("pca_gamma" %in% enrichment_method){
        print("pca_gamma will start.")
        result_pca_gamma=gamma_test(geneset=geneset,value=value)
        result_pca_gamma$two=TRUE
        result_pca_gamma$invert_A=ifelse(result_pca_gamma$r<=0,TRUE,FALSE)
        result_pca_gamma$invert_I=ifelse(result_pca_gamma$r>0,TRUE,FALSE)
        result_pca_gamma$P.Value.A=two2one(result_pca_gamma$P.Value, two = as.logical(result_pca_gamma$two), invert = as.logical(result_pca_gamma$invert_A))
        result_pca_gamma$P.Value.I=two2one(result_pca_gamma$P.Value, two = as.logical(result_pca_gamma$two), invert = as.logical(result_pca_gamma$invert_I))
        result_pca_gamma_A=result_pca_gamma[,c("geneset_id","P.Value.A")]
        result_pca_gamma_I=result_pca_gamma[,c("geneset_id","P.Value.I")]
        colnames(result_pca_gamma_A)[2]=paste0(vector[i],".pca_gamma")
        colnames(result_pca_gamma_I)[2]=paste0(vector[i],".pca_gamma")
        rm(result_pca_gamma)
      }
    }else {
      print("Not using pca.")
    }
    ##########################################aucell
    methods_to_check=c("aucell_pearson","aucell_kendall","aucell_spearman","aucell_lm","aucell_biweight","aucell_distance","aucell_percentage","aucell_blomqvist","aucell_hoeffding","aucell_gamma")
    if (any(methods_to_check %in% enrichment_method)){
      geneset=AUCell_run(exprMat=as.matrix(exp),geneSets=geneSets)
      geneset=geneset@assays@data@listData[["AUC"]]
      geneset=as.data.frame(geneset)
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("aucell_pearson" %in% enrichment_method){
        print("aucell_pearson will start.")
        result_aucell_pearson=pearson_test(geneset=geneset,value=value)
        result_aucell_pearson$two=TRUE
        result_aucell_pearson$invert_A=ifelse(result_aucell_pearson$r<=0,TRUE,FALSE)
        result_aucell_pearson$invert_I=ifelse(result_aucell_pearson$r>0,TRUE,FALSE)
        result_aucell_pearson$P.Value.A=two2one(result_aucell_pearson$P.Value, two = as.logical(result_aucell_pearson$two), invert = as.logical(result_aucell_pearson$invert_A))
        result_aucell_pearson$P.Value.I=two2one(result_aucell_pearson$P.Value, two = as.logical(result_aucell_pearson$two), invert = as.logical(result_aucell_pearson$invert_I))
        result_aucell_pearson_A=result_aucell_pearson[,c("geneset_id","P.Value.A")]
        result_aucell_pearson_I=result_aucell_pearson[,c("geneset_id","P.Value.I")]
        colnames(result_aucell_pearson_A)[2]=paste0(vector[i],".aucell_pearson")
        colnames(result_aucell_pearson_I)[2]=paste0(vector[i],".aucell_pearson")
        rm(result_aucell_pearson)
      }
      if("aucell_kendall" %in% enrichment_method){
        print("aucell_kendall will start.")
        result_aucell_kendall=kendall_test(geneset=geneset,value=value)
        result_aucell_kendall$two=TRUE
        result_aucell_kendall$invert_A=ifelse(result_aucell_kendall$r<=0,TRUE,FALSE)
        result_aucell_kendall$invert_I=ifelse(result_aucell_kendall$r>0,TRUE,FALSE)
        result_aucell_kendall$P.Value.A=two2one(result_aucell_kendall$P.Value, two = as.logical(result_aucell_kendall$two), invert = as.logical(result_aucell_kendall$invert_A))
        result_aucell_kendall$P.Value.I=two2one(result_aucell_kendall$P.Value, two = as.logical(result_aucell_kendall$two), invert = as.logical(result_aucell_kendall$invert_I))
        result_aucell_kendall_A=result_aucell_kendall[,c("geneset_id","P.Value.A")]
        result_aucell_kendall_I=result_aucell_kendall[,c("geneset_id","P.Value.I")]
        colnames(result_aucell_kendall_A)[2]=paste0(vector[i],".aucell_kendall")
        colnames(result_aucell_kendall_I)[2]=paste0(vector[i],".aucell_kendall")
        rm(result_aucell_kendall)
      }
      if("aucell_spearman" %in% enrichment_method){
        print("aucell_spearman will start.")
        result_aucell_spearman=spearman_test(geneset=geneset,value=value)
        result_aucell_spearman$two=TRUE
        result_aucell_spearman$invert_A=ifelse(result_aucell_spearman$r<=0,TRUE,FALSE)
        result_aucell_spearman$invert_I=ifelse(result_aucell_spearman$r>0,TRUE,FALSE)
        result_aucell_spearman$P.Value.A=two2one(result_aucell_spearman$P.Value, two = as.logical(result_aucell_spearman$two), invert = as.logical(result_aucell_spearman$invert_A))
        result_aucell_spearman$P.Value.I=two2one(result_aucell_spearman$P.Value, two = as.logical(result_aucell_spearman$two), invert = as.logical(result_aucell_spearman$invert_I))
        result_aucell_spearman_A=result_aucell_spearman[,c("geneset_id","P.Value.A")]
        result_aucell_spearman_I=result_aucell_spearman[,c("geneset_id","P.Value.I")]
        colnames(result_aucell_spearman_A)[2]=paste0(vector[i],".aucell_spearman")
        colnames(result_aucell_spearman_I)[2]=paste0(vector[i],".aucell_spearman")
        rm(result_aucell_spearman)
      }
      if("aucell_lm" %in% enrichment_method){
        print("aucell_lm will start.")
        result_aucell_lm=lm_test(geneset=geneset,value=value)
        result_aucell_lm$two=TRUE
        result_aucell_lm$invert_A=ifelse(result_aucell_lm$r<=0,TRUE,FALSE)
        result_aucell_lm$invert_I=ifelse(result_aucell_lm$r>0,TRUE,FALSE)
        result_aucell_lm$P.Value.A=two2one(result_aucell_lm$P.Value, two = as.logical(result_aucell_lm$two), invert = as.logical(result_aucell_lm$invert_A))
        result_aucell_lm$P.Value.I=two2one(result_aucell_lm$P.Value, two = as.logical(result_aucell_lm$two), invert = as.logical(result_aucell_lm$invert_I))
        result_aucell_lm_A=result_aucell_lm[,c("geneset_id","P.Value.A")]
        result_aucell_lm_I=result_aucell_lm[,c("geneset_id","P.Value.I")]
        colnames(result_aucell_lm_A)[2]=paste0(vector[i],".aucell_lm")
        colnames(result_aucell_lm_I)[2]=paste0(vector[i],".aucell_lm")
        rm(result_aucell_lm)
      }
      if("aucell_biweight" %in% enrichment_method){
        print("aucell_biweight will start.")
        result_aucell_biweight=biweight_test(geneset=geneset,value=value)
        result_aucell_biweight$two=TRUE
        result_aucell_biweight$invert_A=ifelse(result_aucell_biweight$r<=0,TRUE,FALSE)
        result_aucell_biweight$invert_I=ifelse(result_aucell_biweight$r>0,TRUE,FALSE)
        result_aucell_biweight$P.Value.A=two2one(result_aucell_biweight$P.Value, two = as.logical(result_aucell_biweight$two), invert = as.logical(result_aucell_biweight$invert_A))
        result_aucell_biweight$P.Value.I=two2one(result_aucell_biweight$P.Value, two = as.logical(result_aucell_biweight$two), invert = as.logical(result_aucell_biweight$invert_I))
        result_aucell_biweight_A=result_aucell_biweight[,c("geneset_id","P.Value.A")]
        result_aucell_biweight_I=result_aucell_biweight[,c("geneset_id","P.Value.I")]
        colnames(result_aucell_biweight_A)[2]=paste0(vector[i],".aucell_biweight")
        colnames(result_aucell_biweight_I)[2]=paste0(vector[i],".aucell_biweight")
        rm(result_aucell_biweight)
      }
      if("aucell_distance" %in% enrichment_method){
        print("aucell_distance will start.")
        result_aucell_distance=distance_test(geneset=geneset,value=value)
        result_aucell_distance$two=TRUE
        result_aucell_distance$invert_A=ifelse(result_aucell_distance$r<=0,TRUE,FALSE)
        result_aucell_distance$invert_I=ifelse(result_aucell_distance$r>0,TRUE,FALSE)
        result_aucell_distance$P.Value.A=two2one(result_aucell_distance$P.Value, two = as.logical(result_aucell_distance$two), invert = as.logical(result_aucell_distance$invert_A))
        result_aucell_distance$P.Value.I=two2one(result_aucell_distance$P.Value, two = as.logical(result_aucell_distance$two), invert = as.logical(result_aucell_distance$invert_I))
        result_aucell_distance_A=result_aucell_distance[,c("geneset_id","P.Value.A")]
        result_aucell_distance_I=result_aucell_distance[,c("geneset_id","P.Value.I")]
        colnames(result_aucell_distance_A)[2]=paste0(vector[i],".aucell_distance")
        colnames(result_aucell_distance_I)[2]=paste0(vector[i],".aucell_distance")
        rm(result_aucell_distance)
      }
      if("aucell_percentage" %in% enrichment_method){
        print("aucell_percentage will start.")
        result_aucell_percentage=percentage_test(geneset=geneset,value=value)
        result_aucell_percentage$two=TRUE
        result_aucell_percentage$invert_A=ifelse(result_aucell_percentage$r<=0,TRUE,FALSE)
        result_aucell_percentage$invert_I=ifelse(result_aucell_percentage$r>0,TRUE,FALSE)
        result_aucell_percentage$P.Value.A=two2one(result_aucell_percentage$P.Value, two = as.logical(result_aucell_percentage$two), invert = as.logical(result_aucell_percentage$invert_A))
        result_aucell_percentage$P.Value.I=two2one(result_aucell_percentage$P.Value, two = as.logical(result_aucell_percentage$two), invert = as.logical(result_aucell_percentage$invert_I))
        result_aucell_percentage_A=result_aucell_percentage[,c("geneset_id","P.Value.A")]
        result_aucell_percentage_I=result_aucell_percentage[,c("geneset_id","P.Value.I")]
        colnames(result_aucell_percentage_A)[2]=paste0(vector[i],".aucell_percentage")
        colnames(result_aucell_percentage_I)[2]=paste0(vector[i],".aucell_percentage")
        rm(result_aucell_percentage)
      }
      if("aucell_blomqvist" %in% enrichment_method){
        print("aucell_blomqvist will start.")
        result_aucell_blomqvist=blomqvist_test(geneset=geneset,value=value)
        result_aucell_blomqvist$two=TRUE
        result_aucell_blomqvist$invert_A=ifelse(result_aucell_blomqvist$r<=0,TRUE,FALSE)
        result_aucell_blomqvist$invert_I=ifelse(result_aucell_blomqvist$r>0,TRUE,FALSE)
        result_aucell_blomqvist$P.Value.A=two2one(result_aucell_blomqvist$P.Value, two = as.logical(result_aucell_blomqvist$two), invert = as.logical(result_aucell_blomqvist$invert_A))
        result_aucell_blomqvist$P.Value.I=two2one(result_aucell_blomqvist$P.Value, two = as.logical(result_aucell_blomqvist$two), invert = as.logical(result_aucell_blomqvist$invert_I))
        result_aucell_blomqvist_A=result_aucell_blomqvist[,c("geneset_id","P.Value.A")]
        result_aucell_blomqvist_I=result_aucell_blomqvist[,c("geneset_id","P.Value.I")]
        colnames(result_aucell_blomqvist_A)[2]=paste0(vector[i],".aucell_blomqvist")
        colnames(result_aucell_blomqvist_I)[2]=paste0(vector[i],".aucell_blomqvist")
        rm(result_aucell_blomqvist)
      }
      if("aucell_hoeffding" %in% enrichment_method){
        print("aucell_hoeffding will start.")
        result_aucell_hoeffding=hoeffding_test(geneset=geneset,value=value)
        result_aucell_hoeffding$two=TRUE
        result_aucell_hoeffding$invert_A=ifelse(result_aucell_hoeffding$r<=0,TRUE,FALSE)
        result_aucell_hoeffding$invert_I=ifelse(result_aucell_hoeffding$r>0,TRUE,FALSE)
        result_aucell_hoeffding$P.Value.A=two2one(result_aucell_hoeffding$P.Value, two = as.logical(result_aucell_hoeffding$two), invert = as.logical(result_aucell_hoeffding$invert_A))
        result_aucell_hoeffding$P.Value.I=two2one(result_aucell_hoeffding$P.Value, two = as.logical(result_aucell_hoeffding$two), invert = as.logical(result_aucell_hoeffding$invert_I))
        result_aucell_hoeffding_A=result_aucell_hoeffding[,c("geneset_id","P.Value.A")]
        result_aucell_hoeffding_I=result_aucell_hoeffding[,c("geneset_id","P.Value.I")]
        colnames(result_aucell_hoeffding_A)[2]=paste0(vector[i],".aucell_hoeffding")
        colnames(result_aucell_hoeffding_I)[2]=paste0(vector[i],".aucell_hoeffding")
        rm(result_aucell_hoeffding)
      }
      if("aucell_gamma" %in% enrichment_method){
        print("aucell_gamma will start.")
        result_aucell_gamma=gamma_test(geneset=geneset,value=value)
        result_aucell_gamma$two=TRUE
        result_aucell_gamma$invert_A=ifelse(result_aucell_gamma$r<=0,TRUE,FALSE)
        result_aucell_gamma$invert_I=ifelse(result_aucell_gamma$r>0,TRUE,FALSE)
        result_aucell_gamma$P.Value.A=two2one(result_aucell_gamma$P.Value, two = as.logical(result_aucell_gamma$two), invert = as.logical(result_aucell_gamma$invert_A))
        result_aucell_gamma$P.Value.I=two2one(result_aucell_gamma$P.Value, two = as.logical(result_aucell_gamma$two), invert = as.logical(result_aucell_gamma$invert_I))
        result_aucell_gamma_A=result_aucell_gamma[,c("geneset_id","P.Value.A")]
        result_aucell_gamma_I=result_aucell_gamma[,c("geneset_id","P.Value.I")]
        colnames(result_aucell_gamma_A)[2]=paste0(vector[i],".aucell_gamma")
        colnames(result_aucell_gamma_I)[2]=paste0(vector[i],".aucell_gamma")
        rm(result_aucell_gamma)
      }
    }else {
      print("Not using aucell.")
    }
    ##########################################ucell
    methods_to_check=c("ucell_pearson","ucell_kendall","ucell_spearman","ucell_lm","ucell_biweight","ucell_distance","ucell_percentage","ucell_blomqvist","ucell_hoeffding","ucell_gamma")
    if (any(methods_to_check %in% enrichment_method)){
      geneSets_list=lapply(geneSets, function(x){ x@geneIds })
      names(geneSets_list)=names(geneSets)
      geneset=ScoreSignatures_UCell(as.matrix(exp), geneSets_list)
      geneset=as.data.frame(t(geneset))
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("ucell_pearson" %in% enrichment_method){
        print("ucell_pearson will start.")
        result_ucell_pearson=pearson_test(geneset=geneset,value=value)
        result_ucell_pearson$geneset_id=gsub("_UCell$", "", result_ucell_pearson$geneset_id)
        result_ucell_pearson$two=TRUE
        result_ucell_pearson$invert_A=ifelse(result_ucell_pearson$r<=0,TRUE,FALSE)
        result_ucell_pearson$invert_I=ifelse(result_ucell_pearson$r>0,TRUE,FALSE)
        result_ucell_pearson$P.Value.A=two2one(result_ucell_pearson$P.Value, two = as.logical(result_ucell_pearson$two), invert = as.logical(result_ucell_pearson$invert_A))
        result_ucell_pearson$P.Value.I=two2one(result_ucell_pearson$P.Value, two = as.logical(result_ucell_pearson$two), invert = as.logical(result_ucell_pearson$invert_I))
        result_ucell_pearson_A=result_ucell_pearson[,c("geneset_id","P.Value.A")]
        result_ucell_pearson_I=result_ucell_pearson[,c("geneset_id","P.Value.I")]
        colnames(result_ucell_pearson_A)[2]=paste0(vector[i],".ucell_pearson")
        colnames(result_ucell_pearson_I)[2]=paste0(vector[i],".ucell_pearson")
        rm(result_ucell_pearson)
      }
      if("ucell_kendall" %in% enrichment_method){
        print("ucell_kendall will start.")
        result_ucell_kendall=kendall_test(geneset=geneset,value=value)
        result_ucell_kendall$geneset_id=gsub("_UCell$", "", result_ucell_kendall$geneset_id)
        result_ucell_kendall$two=TRUE
        result_ucell_kendall$invert_A=ifelse(result_ucell_kendall$r<=0,TRUE,FALSE)
        result_ucell_kendall$invert_I=ifelse(result_ucell_kendall$r>0,TRUE,FALSE)
        result_ucell_kendall$P.Value.A=two2one(result_ucell_kendall$P.Value, two = as.logical(result_ucell_kendall$two), invert = as.logical(result_ucell_kendall$invert_A))
        result_ucell_kendall$P.Value.I=two2one(result_ucell_kendall$P.Value, two = as.logical(result_ucell_kendall$two), invert = as.logical(result_ucell_kendall$invert_I))
        result_ucell_kendall_A=result_ucell_kendall[,c("geneset_id","P.Value.A")]
        result_ucell_kendall_I=result_ucell_kendall[,c("geneset_id","P.Value.I")]
        colnames(result_ucell_kendall_A)[2]=paste0(vector[i],".ucell_kendall")
        colnames(result_ucell_kendall_I)[2]=paste0(vector[i],".ucell_kendall")
        rm(result_ucell_kendall)
      }
      if("ucell_spearman" %in% enrichment_method){
        print("ucell_spearman will start.")
        result_ucell_spearman=spearman_test(geneset=geneset,value=value)
        result_ucell_spearman$geneset_id=gsub("_UCell$", "", result_ucell_spearman$geneset_id)
        result_ucell_spearman$two=TRUE
        result_ucell_spearman$invert_A=ifelse(result_ucell_spearman$r<=0,TRUE,FALSE)
        result_ucell_spearman$invert_I=ifelse(result_ucell_spearman$r>0,TRUE,FALSE)
        result_ucell_spearman$P.Value.A=two2one(result_ucell_spearman$P.Value, two = as.logical(result_ucell_spearman$two), invert = as.logical(result_ucell_spearman$invert_A))
        result_ucell_spearman$P.Value.I=two2one(result_ucell_spearman$P.Value, two = as.logical(result_ucell_spearman$two), invert = as.logical(result_ucell_spearman$invert_I))
        result_ucell_spearman_A=result_ucell_spearman[,c("geneset_id","P.Value.A")]
        result_ucell_spearman_I=result_ucell_spearman[,c("geneset_id","P.Value.I")]
        colnames(result_ucell_spearman_A)[2]=paste0(vector[i],".ucell_spearman")
        colnames(result_ucell_spearman_I)[2]=paste0(vector[i],".ucell_spearman")
        rm(result_ucell_spearman)
      }
      if("ucell_lm" %in% enrichment_method){
        print("ucell_lm will start.")
        result_ucell_lm=lm_test(geneset=geneset,value=value)
        result_ucell_lm$geneset_id=gsub("_UCell$", "", result_ucell_lm$geneset_id)
        result_ucell_lm$two=TRUE
        result_ucell_lm$invert_A=ifelse(result_ucell_lm$r<=0,TRUE,FALSE)
        result_ucell_lm$invert_I=ifelse(result_ucell_lm$r>0,TRUE,FALSE)
        result_ucell_lm$P.Value.A=two2one(result_ucell_lm$P.Value, two = as.logical(result_ucell_lm$two), invert = as.logical(result_ucell_lm$invert_A))
        result_ucell_lm$P.Value.I=two2one(result_ucell_lm$P.Value, two = as.logical(result_ucell_lm$two), invert = as.logical(result_ucell_lm$invert_I))
        result_ucell_lm_A=result_ucell_lm[,c("geneset_id","P.Value.A")]
        result_ucell_lm_I=result_ucell_lm[,c("geneset_id","P.Value.I")]
        colnames(result_ucell_lm_A)[2]=paste0(vector[i],".ucell_lm")
        colnames(result_ucell_lm_I)[2]=paste0(vector[i],".ucell_lm")
        rm(result_ucell_lm)
      }
      if("ucell_biweight" %in% enrichment_method){
        print("ucell_biweight will start.")
        result_ucell_biweight=biweight_test(geneset=geneset,value=value)
        result_ucell_biweight$geneset_id=gsub("_UCell$", "", result_ucell_biweight$geneset_id)
        result_ucell_biweight$two=TRUE
        result_ucell_biweight$invert_A=ifelse(result_ucell_biweight$r<=0,TRUE,FALSE)
        result_ucell_biweight$invert_I=ifelse(result_ucell_biweight$r>0,TRUE,FALSE)
        result_ucell_biweight$P.Value.A=two2one(result_ucell_biweight$P.Value, two = as.logical(result_ucell_biweight$two), invert = as.logical(result_ucell_biweight$invert_A))
        result_ucell_biweight$P.Value.I=two2one(result_ucell_biweight$P.Value, two = as.logical(result_ucell_biweight$two), invert = as.logical(result_ucell_biweight$invert_I))
        result_ucell_biweight_A=result_ucell_biweight[,c("geneset_id","P.Value.A")]
        result_ucell_biweight_I=result_ucell_biweight[,c("geneset_id","P.Value.I")]
        colnames(result_ucell_biweight_A)[2]=paste0(vector[i],".ucell_biweight")
        colnames(result_ucell_biweight_I)[2]=paste0(vector[i],".ucell_biweight")
        rm(result_ucell_biweight)
      }
      if("ucell_distance" %in% enrichment_method){
        print("ucell_distance will start.")
        result_ucell_distance=distance_test(geneset=geneset,value=value)
        result_ucell_distance$geneset_id=gsub("_UCell$", "", result_ucell_distance$geneset_id)
        result_ucell_distance$two=TRUE
        result_ucell_distance$invert_A=ifelse(result_ucell_distance$r<=0,TRUE,FALSE)
        result_ucell_distance$invert_I=ifelse(result_ucell_distance$r>0,TRUE,FALSE)
        result_ucell_distance$P.Value.A=two2one(result_ucell_distance$P.Value, two = as.logical(result_ucell_distance$two), invert = as.logical(result_ucell_distance$invert_A))
        result_ucell_distance$P.Value.I=two2one(result_ucell_distance$P.Value, two = as.logical(result_ucell_distance$two), invert = as.logical(result_ucell_distance$invert_I))
        result_ucell_distance_A=result_ucell_distance[,c("geneset_id","P.Value.A")]
        result_ucell_distance_I=result_ucell_distance[,c("geneset_id","P.Value.I")]
        colnames(result_ucell_distance_A)[2]=paste0(vector[i],".ucell_distance")
        colnames(result_ucell_distance_I)[2]=paste0(vector[i],".ucell_distance")
        rm(result_ucell_distance)
      }
      if("ucell_percentage" %in% enrichment_method){
        print("ucell_percentage will start.")
        result_ucell_percentage=percentage_test(geneset=geneset,value=value)
        result_ucell_percentage$geneset_id=gsub("_UCell$", "", result_ucell_percentage$geneset_id)
        result_ucell_percentage$two=TRUE
        result_ucell_percentage$invert_A=ifelse(result_ucell_percentage$r<=0,TRUE,FALSE)
        result_ucell_percentage$invert_I=ifelse(result_ucell_percentage$r>0,TRUE,FALSE)
        result_ucell_percentage$P.Value.A=two2one(result_ucell_percentage$P.Value, two = as.logical(result_ucell_percentage$two), invert = as.logical(result_ucell_percentage$invert_A))
        result_ucell_percentage$P.Value.I=two2one(result_ucell_percentage$P.Value, two = as.logical(result_ucell_percentage$two), invert = as.logical(result_ucell_percentage$invert_I))
        result_ucell_percentage_A=result_ucell_percentage[,c("geneset_id","P.Value.A")]
        result_ucell_percentage_I=result_ucell_percentage[,c("geneset_id","P.Value.I")]
        colnames(result_ucell_percentage_A)[2]=paste0(vector[i],".ucell_percentage")
        colnames(result_ucell_percentage_I)[2]=paste0(vector[i],".ucell_percentage")
        rm(result_ucell_percentage)
      }
      if("ucell_blomqvist" %in% enrichment_method){
        print("ucell_blomqvist will start.")
        result_ucell_blomqvist=blomqvist_test(geneset=geneset,value=value)
        result_ucell_blomqvist$geneset_id=gsub("_UCell$", "", result_ucell_blomqvist$geneset_id)
        result_ucell_blomqvist$two=TRUE
        result_ucell_blomqvist$invert_A=ifelse(result_ucell_blomqvist$r<=0,TRUE,FALSE)
        result_ucell_blomqvist$invert_I=ifelse(result_ucell_blomqvist$r>0,TRUE,FALSE)
        result_ucell_blomqvist$P.Value.A=two2one(result_ucell_blomqvist$P.Value, two = as.logical(result_ucell_blomqvist$two), invert = as.logical(result_ucell_blomqvist$invert_A))
        result_ucell_blomqvist$P.Value.I=two2one(result_ucell_blomqvist$P.Value, two = as.logical(result_ucell_blomqvist$two), invert = as.logical(result_ucell_blomqvist$invert_I))
        result_ucell_blomqvist_A=result_ucell_blomqvist[,c("geneset_id","P.Value.A")]
        result_ucell_blomqvist_I=result_ucell_blomqvist[,c("geneset_id","P.Value.I")]
        colnames(result_ucell_blomqvist_A)[2]=paste0(vector[i],".ucell_blomqvist")
        colnames(result_ucell_blomqvist_I)[2]=paste0(vector[i],".ucell_blomqvist")
        rm(result_ucell_blomqvist)
      }
      if("ucell_hoeffding" %in% enrichment_method){
        print("ucell_hoeffding will start.")
        result_ucell_hoeffding=hoeffding_test(geneset=geneset,value=value)
        result_ucell_hoeffding$geneset_id=gsub("_UCell$", "", result_ucell_hoeffding$geneset_id)
        result_ucell_hoeffding$two=TRUE
        result_ucell_hoeffding$invert_A=ifelse(result_ucell_hoeffding$r<=0,TRUE,FALSE)
        result_ucell_hoeffding$invert_I=ifelse(result_ucell_hoeffding$r>0,TRUE,FALSE)
        result_ucell_hoeffding$P.Value.A=two2one(result_ucell_hoeffding$P.Value, two = as.logical(result_ucell_hoeffding$two), invert = as.logical(result_ucell_hoeffding$invert_A))
        result_ucell_hoeffding$P.Value.I=two2one(result_ucell_hoeffding$P.Value, two = as.logical(result_ucell_hoeffding$two), invert = as.logical(result_ucell_hoeffding$invert_I))
        result_ucell_hoeffding_A=result_ucell_hoeffding[,c("geneset_id","P.Value.A")]
        result_ucell_hoeffding_I=result_ucell_hoeffding[,c("geneset_id","P.Value.I")]
        colnames(result_ucell_hoeffding_A)[2]=paste0(vector[i],".ucell_hoeffding")
        colnames(result_ucell_hoeffding_I)[2]=paste0(vector[i],".ucell_hoeffding")
        rm(result_ucell_hoeffding)
      }
      if("ucell_gamma" %in% enrichment_method){
        print("ucell_gamma will start.")
        result_ucell_gamma=gamma_test(geneset=geneset,value=value)
        result_ucell_gamma$geneset_id=gsub("_UCell$", "", result_ucell_gamma$geneset_id)
        result_ucell_gamma$two=TRUE
        result_ucell_gamma$invert_A=ifelse(result_ucell_gamma$r<=0,TRUE,FALSE)
        result_ucell_gamma$invert_I=ifelse(result_ucell_gamma$r>0,TRUE,FALSE)
        result_ucell_gamma$P.Value.A=two2one(result_ucell_gamma$P.Value, two = as.logical(result_ucell_gamma$two), invert = as.logical(result_ucell_gamma$invert_A))
        result_ucell_gamma$P.Value.I=two2one(result_ucell_gamma$P.Value, two = as.logical(result_ucell_gamma$two), invert = as.logical(result_ucell_gamma$invert_I))
        result_ucell_gamma_A=result_ucell_gamma[,c("geneset_id","P.Value.A")]
        result_ucell_gamma_I=result_ucell_gamma[,c("geneset_id","P.Value.I")]
        colnames(result_ucell_gamma_A)[2]=paste0(vector[i],".ucell_gamma")
        colnames(result_ucell_gamma_I)[2]=paste0(vector[i],".ucell_gamma")
        rm(result_ucell_gamma)
      }
    }else {
      print("Not using ucell.")
    }
    ##########################################singscore
    methods_to_check=c("singscore_pearson","singscore_kendall","singscore_spearman","singscore_lm","singscore_biweight","singscore_distance","singscore_percentage","singscore_blomqvist","singscore_hoeffding","singscore_gamma")
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
      if("singscore_pearson" %in% enrichment_method){
        print("singscore_pearson will start.")
        result_singscore_pearson=pearson_test(geneset=geneset,value=value)
        result_singscore_pearson$two=TRUE
        result_singscore_pearson$invert_A=ifelse(result_singscore_pearson$r<=0,TRUE,FALSE)
        result_singscore_pearson$invert_I=ifelse(result_singscore_pearson$r>0,TRUE,FALSE)
        result_singscore_pearson$P.Value.A=two2one(result_singscore_pearson$P.Value, two = as.logical(result_singscore_pearson$two), invert = as.logical(result_singscore_pearson$invert_A))
        result_singscore_pearson$P.Value.I=two2one(result_singscore_pearson$P.Value, two = as.logical(result_singscore_pearson$two), invert = as.logical(result_singscore_pearson$invert_I))
        result_singscore_pearson_A=result_singscore_pearson[,c("geneset_id","P.Value.A")]
        result_singscore_pearson_I=result_singscore_pearson[,c("geneset_id","P.Value.I")]
        colnames(result_singscore_pearson_A)[2]=paste0(vector[i],".singscore_pearson")
        colnames(result_singscore_pearson_I)[2]=paste0(vector[i],".singscore_pearson")
        rm(result_singscore_pearson)
      }
      if("singscore_kendall" %in% enrichment_method){
        print("singscore_kendall will start.")
        result_singscore_kendall=kendall_test(geneset=geneset,value=value)
        result_singscore_kendall$two=TRUE
        result_singscore_kendall$invert_A=ifelse(result_singscore_kendall$r<=0,TRUE,FALSE)
        result_singscore_kendall$invert_I=ifelse(result_singscore_kendall$r>0,TRUE,FALSE)
        result_singscore_kendall$P.Value.A=two2one(result_singscore_kendall$P.Value, two = as.logical(result_singscore_kendall$two), invert = as.logical(result_singscore_kendall$invert_A))
        result_singscore_kendall$P.Value.I=two2one(result_singscore_kendall$P.Value, two = as.logical(result_singscore_kendall$two), invert = as.logical(result_singscore_kendall$invert_I))
        result_singscore_kendall_A=result_singscore_kendall[,c("geneset_id","P.Value.A")]
        result_singscore_kendall_I=result_singscore_kendall[,c("geneset_id","P.Value.I")]
        colnames(result_singscore_kendall_A)[2]=paste0(vector[i],".singscore_kendall")
        colnames(result_singscore_kendall_I)[2]=paste0(vector[i],".singscore_kendall")
        rm(result_singscore_kendall)
      }
      if("singscore_spearman" %in% enrichment_method){
        print("singscore_spearman will start.")
        result_singscore_spearman=spearman_test(geneset=geneset,value=value)
        result_singscore_spearman$two=TRUE
        result_singscore_spearman$invert_A=ifelse(result_singscore_spearman$r<=0,TRUE,FALSE)
        result_singscore_spearman$invert_I=ifelse(result_singscore_spearman$r>0,TRUE,FALSE)
        result_singscore_spearman$P.Value.A=two2one(result_singscore_spearman$P.Value, two = as.logical(result_singscore_spearman$two), invert = as.logical(result_singscore_spearman$invert_A))
        result_singscore_spearman$P.Value.I=two2one(result_singscore_spearman$P.Value, two = as.logical(result_singscore_spearman$two), invert = as.logical(result_singscore_spearman$invert_I))
        result_singscore_spearman_A=result_singscore_spearman[,c("geneset_id","P.Value.A")]
        result_singscore_spearman_I=result_singscore_spearman[,c("geneset_id","P.Value.I")]
        colnames(result_singscore_spearman_A)[2]=paste0(vector[i],".singscore_spearman")
        colnames(result_singscore_spearman_I)[2]=paste0(vector[i],".singscore_spearman")
        rm(result_singscore_spearman)
      }
      if("singscore_lm" %in% enrichment_method){
        print("singscore_lm will start.")
        result_singscore_lm=lm_test(geneset=geneset,value=value)
        result_singscore_lm$two=TRUE
        result_singscore_lm$invert_A=ifelse(result_singscore_lm$r<=0,TRUE,FALSE)
        result_singscore_lm$invert_I=ifelse(result_singscore_lm$r>0,TRUE,FALSE)
        result_singscore_lm$P.Value.A=two2one(result_singscore_lm$P.Value, two = as.logical(result_singscore_lm$two), invert = as.logical(result_singscore_lm$invert_A))
        result_singscore_lm$P.Value.I=two2one(result_singscore_lm$P.Value, two = as.logical(result_singscore_lm$two), invert = as.logical(result_singscore_lm$invert_I))
        result_singscore_lm_A=result_singscore_lm[,c("geneset_id","P.Value.A")]
        result_singscore_lm_I=result_singscore_lm[,c("geneset_id","P.Value.I")]
        colnames(result_singscore_lm_A)[2]=paste0(vector[i],".singscore_lm")
        colnames(result_singscore_lm_I)[2]=paste0(vector[i],".singscore_lm")
        rm(result_singscore_lm)
      }
      if("singscore_biweight" %in% enrichment_method){
        print("singscore_biweight will start.")
        result_singscore_biweight=biweight_test(geneset=geneset,value=value)
        result_singscore_biweight$two=TRUE
        result_singscore_biweight$invert_A=ifelse(result_singscore_biweight$r<=0,TRUE,FALSE)
        result_singscore_biweight$invert_I=ifelse(result_singscore_biweight$r>0,TRUE,FALSE)
        result_singscore_biweight$P.Value.A=two2one(result_singscore_biweight$P.Value, two = as.logical(result_singscore_biweight$two), invert = as.logical(result_singscore_biweight$invert_A))
        result_singscore_biweight$P.Value.I=two2one(result_singscore_biweight$P.Value, two = as.logical(result_singscore_biweight$two), invert = as.logical(result_singscore_biweight$invert_I))
        result_singscore_biweight_A=result_singscore_biweight[,c("geneset_id","P.Value.A")]
        result_singscore_biweight_I=result_singscore_biweight[,c("geneset_id","P.Value.I")]
        colnames(result_singscore_biweight_A)[2]=paste0(vector[i],".singscore_biweight")
        colnames(result_singscore_biweight_I)[2]=paste0(vector[i],".singscore_biweight")
        rm(result_singscore_biweight)
      }
      if("singscore_distance" %in% enrichment_method){
        print("singscore_distance will start.")
        result_singscore_distance=distance_test(geneset=geneset,value=value)
        result_singscore_distance$two=TRUE
        result_singscore_distance$invert_A=ifelse(result_singscore_distance$r<=0,TRUE,FALSE)
        result_singscore_distance$invert_I=ifelse(result_singscore_distance$r>0,TRUE,FALSE)
        result_singscore_distance$P.Value.A=two2one(result_singscore_distance$P.Value, two = as.logical(result_singscore_distance$two), invert = as.logical(result_singscore_distance$invert_A))
        result_singscore_distance$P.Value.I=two2one(result_singscore_distance$P.Value, two = as.logical(result_singscore_distance$two), invert = as.logical(result_singscore_distance$invert_I))
        result_singscore_distance_A=result_singscore_distance[,c("geneset_id","P.Value.A")]
        result_singscore_distance_I=result_singscore_distance[,c("geneset_id","P.Value.I")]
        colnames(result_singscore_distance_A)[2]=paste0(vector[i],".singscore_distance")
        colnames(result_singscore_distance_I)[2]=paste0(vector[i],".singscore_distance")
        rm(result_singscore_distance)
      }
      if("singscore_percentage" %in% enrichment_method){
        print("singscore_percentage will start.")
        result_singscore_percentage=percentage_test(geneset=geneset,value=value)
        result_singscore_percentage$two=TRUE
        result_singscore_percentage$invert_A=ifelse(result_singscore_percentage$r<=0,TRUE,FALSE)
        result_singscore_percentage$invert_I=ifelse(result_singscore_percentage$r>0,TRUE,FALSE)
        result_singscore_percentage$P.Value.A=two2one(result_singscore_percentage$P.Value, two = as.logical(result_singscore_percentage$two), invert = as.logical(result_singscore_percentage$invert_A))
        result_singscore_percentage$P.Value.I=two2one(result_singscore_percentage$P.Value, two = as.logical(result_singscore_percentage$two), invert = as.logical(result_singscore_percentage$invert_I))
        result_singscore_percentage_A=result_singscore_percentage[,c("geneset_id","P.Value.A")]
        result_singscore_percentage_I=result_singscore_percentage[,c("geneset_id","P.Value.I")]
        colnames(result_singscore_percentage_A)[2]=paste0(vector[i],".singscore_percentage")
        colnames(result_singscore_percentage_I)[2]=paste0(vector[i],".singscore_percentage")
        rm(result_singscore_percentage)
      }
      if("singscore_blomqvist" %in% enrichment_method){
        print("singscore_blomqvist will start.")
        result_singscore_blomqvist=blomqvist_test(geneset=geneset,value=value)
        result_singscore_blomqvist$two=TRUE
        result_singscore_blomqvist$invert_A=ifelse(result_singscore_blomqvist$r<=0,TRUE,FALSE)
        result_singscore_blomqvist$invert_I=ifelse(result_singscore_blomqvist$r>0,TRUE,FALSE)
        result_singscore_blomqvist$P.Value.A=two2one(result_singscore_blomqvist$P.Value, two = as.logical(result_singscore_blomqvist$two), invert = as.logical(result_singscore_blomqvist$invert_A))
        result_singscore_blomqvist$P.Value.I=two2one(result_singscore_blomqvist$P.Value, two = as.logical(result_singscore_blomqvist$two), invert = as.logical(result_singscore_blomqvist$invert_I))
        result_singscore_blomqvist_A=result_singscore_blomqvist[,c("geneset_id","P.Value.A")]
        result_singscore_blomqvist_I=result_singscore_blomqvist[,c("geneset_id","P.Value.I")]
        colnames(result_singscore_blomqvist_A)[2]=paste0(vector[i],".singscore_blomqvist")
        colnames(result_singscore_blomqvist_I)[2]=paste0(vector[i],".singscore_blomqvist")
        rm(result_singscore_blomqvist)
      }
      if("singscore_hoeffding" %in% enrichment_method){
        print("singscore_hoeffding will start.")
        result_singscore_hoeffding=hoeffding_test(geneset=geneset,value=value)
        result_singscore_hoeffding$two=TRUE
        result_singscore_hoeffding$invert_A=ifelse(result_singscore_hoeffding$r<=0,TRUE,FALSE)
        result_singscore_hoeffding$invert_I=ifelse(result_singscore_hoeffding$r>0,TRUE,FALSE)
        result_singscore_hoeffding$P.Value.A=two2one(result_singscore_hoeffding$P.Value, two = as.logical(result_singscore_hoeffding$two), invert = as.logical(result_singscore_hoeffding$invert_A))
        result_singscore_hoeffding$P.Value.I=two2one(result_singscore_hoeffding$P.Value, two = as.logical(result_singscore_hoeffding$two), invert = as.logical(result_singscore_hoeffding$invert_I))
        result_singscore_hoeffding_A=result_singscore_hoeffding[,c("geneset_id","P.Value.A")]
        result_singscore_hoeffding_I=result_singscore_hoeffding[,c("geneset_id","P.Value.I")]
        colnames(result_singscore_hoeffding_A)[2]=paste0(vector[i],".singscore_hoeffding")
        colnames(result_singscore_hoeffding_I)[2]=paste0(vector[i],".singscore_hoeffding")
        rm(result_singscore_hoeffding)
      }
      if("singscore_gamma" %in% enrichment_method){
        print("singscore_gamma will start.")
        result_singscore_gamma=gamma_test(geneset=geneset,value=value)
        result_singscore_gamma$two=TRUE
        result_singscore_gamma$invert_A=ifelse(result_singscore_gamma$r<=0,TRUE,FALSE)
        result_singscore_gamma$invert_I=ifelse(result_singscore_gamma$r>0,TRUE,FALSE)
        result_singscore_gamma$P.Value.A=two2one(result_singscore_gamma$P.Value, two = as.logical(result_singscore_gamma$two), invert = as.logical(result_singscore_gamma$invert_A))
        result_singscore_gamma$P.Value.I=two2one(result_singscore_gamma$P.Value, two = as.logical(result_singscore_gamma$two), invert = as.logical(result_singscore_gamma$invert_I))
        result_singscore_gamma_A=result_singscore_gamma[,c("geneset_id","P.Value.A")]
        result_singscore_gamma_I=result_singscore_gamma[,c("geneset_id","P.Value.I")]
        colnames(result_singscore_gamma_A)[2]=paste0(vector[i],".singscore_gamma")
        colnames(result_singscore_gamma_I)[2]=paste0(vector[i],".singscore_gamma")
        rm(result_singscore_gamma)
      }
    }else {
      print("Not using singscore.")
    }
    ##########################################median
    methods_to_check=c("median_pearson","median_kendall","median_spearman","median_lm","median_biweight","median_distance","median_percentage","median_blomqvist","median_hoeffding","median_gamma")
    if (any(methods_to_check %in% enrichment_method)){
      geneset=Median_evaluation(exp=exp,geneSets=geneSets)
      geneset=geneset[apply(geneset, 1, var) != 0, ]
      if("median_pearson" %in% enrichment_method){
        print("median_pearson will start.")
        result_median_pearson=pearson_test(geneset=geneset,value=value)
        result_median_pearson$two=TRUE
        result_median_pearson$invert_A=ifelse(result_median_pearson$r<=0,TRUE,FALSE)
        result_median_pearson$invert_I=ifelse(result_median_pearson$r>0,TRUE,FALSE)
        result_median_pearson$P.Value.A=two2one(result_median_pearson$P.Value, two = as.logical(result_median_pearson$two), invert = as.logical(result_median_pearson$invert_A))
        result_median_pearson$P.Value.I=two2one(result_median_pearson$P.Value, two = as.logical(result_median_pearson$two), invert = as.logical(result_median_pearson$invert_I))
        result_median_pearson_A=result_median_pearson[,c("geneset_id","P.Value.A")]
        result_median_pearson_I=result_median_pearson[,c("geneset_id","P.Value.I")]
        colnames(result_median_pearson_A)[2]=paste0(vector[i],".median_pearson")
        colnames(result_median_pearson_I)[2]=paste0(vector[i],".median_pearson")
        rm(result_median_pearson)
      }
      if("median_kendall" %in% enrichment_method){
        print("median_kendall will start.")
        result_median_kendall=kendall_test(geneset=geneset,value=value)
        result_median_kendall$two=TRUE
        result_median_kendall$invert_A=ifelse(result_median_kendall$r<=0,TRUE,FALSE)
        result_median_kendall$invert_I=ifelse(result_median_kendall$r>0,TRUE,FALSE)
        result_median_kendall$P.Value.A=two2one(result_median_kendall$P.Value, two = as.logical(result_median_kendall$two), invert = as.logical(result_median_kendall$invert_A))
        result_median_kendall$P.Value.I=two2one(result_median_kendall$P.Value, two = as.logical(result_median_kendall$two), invert = as.logical(result_median_kendall$invert_I))
        result_median_kendall_A=result_median_kendall[,c("geneset_id","P.Value.A")]
        result_median_kendall_I=result_median_kendall[,c("geneset_id","P.Value.I")]
        colnames(result_median_kendall_A)[2]=paste0(vector[i],".median_kendall")
        colnames(result_median_kendall_I)[2]=paste0(vector[i],".median_kendall")
        rm(result_median_kendall)
      }
      if("median_spearman" %in% enrichment_method){
        print("median_spearman will start.")
        result_median_spearman=spearman_test(geneset=geneset,value=value)
        result_median_spearman$two=TRUE
        result_median_spearman$invert_A=ifelse(result_median_spearman$r<=0,TRUE,FALSE)
        result_median_spearman$invert_I=ifelse(result_median_spearman$r>0,TRUE,FALSE)
        result_median_spearman$P.Value.A=two2one(result_median_spearman$P.Value, two = as.logical(result_median_spearman$two), invert = as.logical(result_median_spearman$invert_A))
        result_median_spearman$P.Value.I=two2one(result_median_spearman$P.Value, two = as.logical(result_median_spearman$two), invert = as.logical(result_median_spearman$invert_I))
        result_median_spearman_A=result_median_spearman[,c("geneset_id","P.Value.A")]
        result_median_spearman_I=result_median_spearman[,c("geneset_id","P.Value.I")]
        colnames(result_median_spearman_A)[2]=paste0(vector[i],".median_spearman")
        colnames(result_median_spearman_I)[2]=paste0(vector[i],".median_spearman")
        rm(result_median_spearman)
      }
      if("median_lm" %in% enrichment_method){
        print("median_lm will start.")
        result_median_lm=lm_test(geneset=geneset,value=value)
        result_median_lm$two=TRUE
        result_median_lm$invert_A=ifelse(result_median_lm$r<=0,TRUE,FALSE)
        result_median_lm$invert_I=ifelse(result_median_lm$r>0,TRUE,FALSE)
        result_median_lm$P.Value.A=two2one(result_median_lm$P.Value, two = as.logical(result_median_lm$two), invert = as.logical(result_median_lm$invert_A))
        result_median_lm$P.Value.I=two2one(result_median_lm$P.Value, two = as.logical(result_median_lm$two), invert = as.logical(result_median_lm$invert_I))
        result_median_lm_A=result_median_lm[,c("geneset_id","P.Value.A")]
        result_median_lm_I=result_median_lm[,c("geneset_id","P.Value.I")]
        colnames(result_median_lm_A)[2]=paste0(vector[i],".median_lm")
        colnames(result_median_lm_I)[2]=paste0(vector[i],".median_lm")
        rm(result_median_lm)
      }
      if("median_biweight" %in% enrichment_method){
        print("median_biweight will start.")
        result_median_biweight=biweight_test(geneset=geneset,value=value)
        result_median_biweight$two=TRUE
        result_median_biweight$invert_A=ifelse(result_median_biweight$r<=0,TRUE,FALSE)
        result_median_biweight$invert_I=ifelse(result_median_biweight$r>0,TRUE,FALSE)
        result_median_biweight$P.Value.A=two2one(result_median_biweight$P.Value, two = as.logical(result_median_biweight$two), invert = as.logical(result_median_biweight$invert_A))
        result_median_biweight$P.Value.I=two2one(result_median_biweight$P.Value, two = as.logical(result_median_biweight$two), invert = as.logical(result_median_biweight$invert_I))
        result_median_biweight_A=result_median_biweight[,c("geneset_id","P.Value.A")]
        result_median_biweight_I=result_median_biweight[,c("geneset_id","P.Value.I")]
        colnames(result_median_biweight_A)[2]=paste0(vector[i],".median_biweight")
        colnames(result_median_biweight_I)[2]=paste0(vector[i],".median_biweight")
        rm(result_median_biweight)
      }
      if("median_distance" %in% enrichment_method){
        print("median_distance will start.")
        result_median_distance=distance_test(geneset=geneset,value=value)
        result_median_distance$two=TRUE
        result_median_distance$invert_A=ifelse(result_median_distance$r<=0,TRUE,FALSE)
        result_median_distance$invert_I=ifelse(result_median_distance$r>0,TRUE,FALSE)
        result_median_distance$P.Value.A=two2one(result_median_distance$P.Value, two = as.logical(result_median_distance$two), invert = as.logical(result_median_distance$invert_A))
        result_median_distance$P.Value.I=two2one(result_median_distance$P.Value, two = as.logical(result_median_distance$two), invert = as.logical(result_median_distance$invert_I))
        result_median_distance_A=result_median_distance[,c("geneset_id","P.Value.A")]
        result_median_distance_I=result_median_distance[,c("geneset_id","P.Value.I")]
        colnames(result_median_distance_A)[2]=paste0(vector[i],".median_distance")
        colnames(result_median_distance_I)[2]=paste0(vector[i],".median_distance")
        rm(result_median_distance)
      }
      if("median_percentage" %in% enrichment_method){
        print("median_percentage will start.")
        result_median_percentage=percentage_test(geneset=geneset,value=value)
        result_median_percentage$two=TRUE
        result_median_percentage$invert_A=ifelse(result_median_percentage$r<=0,TRUE,FALSE)
        result_median_percentage$invert_I=ifelse(result_median_percentage$r>0,TRUE,FALSE)
        result_median_percentage$P.Value.A=two2one(result_median_percentage$P.Value, two = as.logical(result_median_percentage$two), invert = as.logical(result_median_percentage$invert_A))
        result_median_percentage$P.Value.I=two2one(result_median_percentage$P.Value, two = as.logical(result_median_percentage$two), invert = as.logical(result_median_percentage$invert_I))
        result_median_percentage_A=result_median_percentage[,c("geneset_id","P.Value.A")]
        result_median_percentage_I=result_median_percentage[,c("geneset_id","P.Value.I")]
        colnames(result_median_percentage_A)[2]=paste0(vector[i],".median_percentage")
        colnames(result_median_percentage_I)[2]=paste0(vector[i],".median_percentage")
        rm(result_median_percentage)
      }
      if("median_blomqvist" %in% enrichment_method){
        print("median_blomqvist will start.")
        result_median_blomqvist=blomqvist_test(geneset=geneset,value=value)
        result_median_blomqvist$two=TRUE
        result_median_blomqvist$invert_A=ifelse(result_median_blomqvist$r<=0,TRUE,FALSE)
        result_median_blomqvist$invert_I=ifelse(result_median_blomqvist$r>0,TRUE,FALSE)
        result_median_blomqvist$P.Value.A=two2one(result_median_blomqvist$P.Value, two = as.logical(result_median_blomqvist$two), invert = as.logical(result_median_blomqvist$invert_A))
        result_median_blomqvist$P.Value.I=two2one(result_median_blomqvist$P.Value, two = as.logical(result_median_blomqvist$two), invert = as.logical(result_median_blomqvist$invert_I))
        result_median_blomqvist_A=result_median_blomqvist[,c("geneset_id","P.Value.A")]
        result_median_blomqvist_I=result_median_blomqvist[,c("geneset_id","P.Value.I")]
        colnames(result_median_blomqvist_A)[2]=paste0(vector[i],".median_blomqvist")
        colnames(result_median_blomqvist_I)[2]=paste0(vector[i],".median_blomqvist")
        rm(result_median_blomqvist)
      }
      if("median_hoeffding" %in% enrichment_method){
        print("median_hoeffding will start.")
        result_median_hoeffding=hoeffding_test(geneset=geneset,value=value)
        result_median_hoeffding$two=TRUE
        result_median_hoeffding$invert_A=ifelse(result_median_hoeffding$r<=0,TRUE,FALSE)
        result_median_hoeffding$invert_I=ifelse(result_median_hoeffding$r>0,TRUE,FALSE)
        result_median_hoeffding$P.Value.A=two2one(result_median_hoeffding$P.Value, two = as.logical(result_median_hoeffding$two), invert = as.logical(result_median_hoeffding$invert_A))
        result_median_hoeffding$P.Value.I=two2one(result_median_hoeffding$P.Value, two = as.logical(result_median_hoeffding$two), invert = as.logical(result_median_hoeffding$invert_I))
        result_median_hoeffding_A=result_median_hoeffding[,c("geneset_id","P.Value.A")]
        result_median_hoeffding_I=result_median_hoeffding[,c("geneset_id","P.Value.I")]
        colnames(result_median_hoeffding_A)[2]=paste0(vector[i],".median_hoeffding")
        colnames(result_median_hoeffding_I)[2]=paste0(vector[i],".median_hoeffding")
        rm(result_median_hoeffding)
      }
      if("median_gamma" %in% enrichment_method){
        print("median_gamma will start.")
        result_median_gamma=gamma_test(geneset=geneset,value=value)
        result_median_gamma$two=TRUE
        result_median_gamma$invert_A=ifelse(result_median_gamma$r<=0,TRUE,FALSE)
        result_median_gamma$invert_I=ifelse(result_median_gamma$r>0,TRUE,FALSE)
        result_median_gamma$P.Value.A=two2one(result_median_gamma$P.Value, two = as.logical(result_median_gamma$two), invert = as.logical(result_median_gamma$invert_A))
        result_median_gamma$P.Value.I=two2one(result_median_gamma$P.Value, two = as.logical(result_median_gamma$two), invert = as.logical(result_median_gamma$invert_I))
        result_median_gamma_A=result_median_gamma[,c("geneset_id","P.Value.A")]
        result_median_gamma_I=result_median_gamma[,c("geneset_id","P.Value.I")]
        colnames(result_median_gamma_A)[2]=paste0(vector[i],".median_gamma")
        colnames(result_median_gamma_I)[2]=paste0(vector[i],".median_gamma")
        rm(result_median_gamma)
      }
    }else {
      print("Not using median.")
    }
    ##########################################fgsea/ora
    methods_to_check=c("pearson_fgsea","kendall_fgsea","spearman_fgsea","lm_fgsea","biweight_fgsea","distance_fgsea","percentage_fgsea","blomqvist_fgsea","hoeffding_fgsea","gamma_fgsea",
                       "pearson_ora","kendall_ora","spearman_ora","lm_ora","biweight_ora","distance_ora","percentage_ora","blomqvist_ora","hoeffding_ora","gamma_ora")
    if (any(methods_to_check %in% enrichment_method)){
      if (any(c("pearson_fgsea","pearson_ora") %in% enrichment_method)){
        exp_fGSEA=pearson_test(geneset=exp,value=value)
        if("pearson_fgsea"%in% enrichment_method){
          print("pearson_fgsea will start.")
          alldiff=exp_fGSEA[,c("geneset_id","P.Value","r")]
          min_non_zero=min(alldiff$P.Value[alldiff$P.Value > 0], na.rm = TRUE)
          alldiff$P.Value[alldiff$P.Value == 0]=min_non_zero * 0.1
          alldiff$id=-log10(alldiff$P.Value)*alldiff$r
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
          result_pearson_fgsea=result[,c("pathway","pval","ES")]
          colnames(result_pearson_fgsea)=c("geneset_id","P.Value","ES")
          result_pearson_fgsea$geneset_id=gsub("\\.", "'", result_pearson_fgsea$geneset_id)
          result_pearson_fgsea=na.omit(result_pearson_fgsea)
          result_pearson_fgsea$two=TRUE
          result_pearson_fgsea$invert_A=ifelse(result_pearson_fgsea$ES<0,TRUE,FALSE)
          result_pearson_fgsea$invert_I=ifelse(result_pearson_fgsea$ES>0,TRUE,FALSE)
          result_pearson_fgsea$P.Value.A=two2one(result_pearson_fgsea$P.Value, two = as.logical(result_pearson_fgsea$two), invert = as.logical(result_pearson_fgsea$invert_A))
          result_pearson_fgsea$P.Value.I=two2one(result_pearson_fgsea$P.Value, two = as.logical(result_pearson_fgsea$two), invert = as.logical(result_pearson_fgsea$invert_I))
          result_pearson_fgsea_A=result_pearson_fgsea[,c("geneset_id","P.Value.A")]
          result_pearson_fgsea_I=result_pearson_fgsea[,c("geneset_id","P.Value.I")]
          colnames(result_pearson_fgsea_A)[2]=paste0(vector[i],".pearson_fgsea")
          colnames(result_pearson_fgsea_I)[2]=paste0(vector[i],".pearson_fgsea")
          rm(result_pearson_fgsea)
        }
        if("pearson_ora"%in% enrichment_method){
          print("pearson_ora will start.")
          result_UP=exp_fGSEA[exp_fGSEA$r>=0,]
          result_DN=exp_fGSEA[exp_fGSEA$r<0,]
          result_UP=result_UP[order(result_UP$P.Value, -result_UP$r), ]
          result_DN=result_DN[order(result_DN$P.Value, result_DN$r), ]
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
            result_pearson_ora_A=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_pearson_ora_A=enrich_H@result
            result_pearson_ora_A=result_pearson_ora_A[,c("Description","pvalue")]
          }
          colnames(result_pearson_ora_A)=c("geneset_id","P.Value_H")
          result_pearson_ora_A$geneset_id=gsub("\\.", "'", result_pearson_ora_A$geneset_id)
          enrich_L=enricher(gene_L, TERM2GENE = canonical_pathways,minGSSize=min.sz,maxGSSize=max.sz,
                            pvalueCutoff = 1,qvalueCutoff = 1)
          if (is.null(enrich_L)) {
            result_pearson_ora_I=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_pearson_ora_I=enrich_L@result
            result_pearson_ora_I=result_pearson_ora_I[,c("Description","pvalue")]
          }
          colnames(result_pearson_ora_I)=c("geneset_id","P.Value_L")
          result_pearson_ora_I$geneset_id=gsub("\\.", "'", result_pearson_ora_I$geneset_id)
          colnames(result_pearson_ora_A)[2]=paste0(vector[i],".pearson_ora")
          colnames(result_pearson_ora_I)[2]=paste0(vector[i],".pearson_ora")
        }
      }
      if (any(c("kendall_fgsea","kendall_ora") %in% enrichment_method)){
        exp_fGSEA=kendall_test(geneset=exp,value=value)
        if("kendall_fgsea" %in% enrichment_method){
          print("kendall_fgsea will start.")
          alldiff=exp_fGSEA[,c("geneset_id","P.Value","r")]
          min_non_zero=min(alldiff$P.Value[alldiff$P.Value > 0], na.rm = TRUE)
          alldiff$P.Value[alldiff$P.Value == 0]=min_non_zero * 0.1
          alldiff$id=-log10(alldiff$P.Value)*alldiff$r
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
          result_kendall_fgsea=result[,c("pathway","pval","ES")]
          colnames(result_kendall_fgsea)=c("geneset_id","P.Value","ES")
          result_kendall_fgsea$geneset_id=gsub("\\.", "'", result_kendall_fgsea$geneset_id)
          result_kendall_fgsea=na.omit(result_kendall_fgsea)
          result_kendall_fgsea$two=TRUE
          result_kendall_fgsea$invert_A=ifelse(result_kendall_fgsea$ES<0,TRUE,FALSE)
          result_kendall_fgsea$invert_I=ifelse(result_kendall_fgsea$ES>0,TRUE,FALSE)
          result_kendall_fgsea$P.Value.A=two2one(result_kendall_fgsea$P.Value, two = as.logical(result_kendall_fgsea$two), invert = as.logical(result_kendall_fgsea$invert_A))
          result_kendall_fgsea$P.Value.I=two2one(result_kendall_fgsea$P.Value, two = as.logical(result_kendall_fgsea$two), invert = as.logical(result_kendall_fgsea$invert_I))
          result_kendall_fgsea_A=result_kendall_fgsea[,c("geneset_id","P.Value.A")]
          result_kendall_fgsea_I=result_kendall_fgsea[,c("geneset_id","P.Value.I")]
          colnames(result_kendall_fgsea_A)[2]=paste0(vector[i],".kendall_fgsea")
          colnames(result_kendall_fgsea_I)[2]=paste0(vector[i],".kendall_fgsea")
          rm(result_kendall_fgsea)
        }
        if("kendall_ora" %in% enrichment_method){
          print("kendall_ora will start.")
          result_UP=exp_fGSEA[exp_fGSEA$r>=0,]
          result_DN=exp_fGSEA[exp_fGSEA$r<0,]
          result_UP=result_UP[order(result_UP$P.Value, -result_UP$r), ]
          result_DN=result_DN[order(result_DN$P.Value, result_DN$r), ]
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
            result_kendall_ora_A=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_kendall_ora_A=enrich_H@result
            result_kendall_ora_A=result_kendall_ora_A[,c("Description","pvalue")]
          }
          colnames(result_kendall_ora_A)=c("geneset_id","P.Value_H")
          result_kendall_ora_A$geneset_id=gsub("\\.", "'", result_kendall_ora_A$geneset_id)
          enrich_L=enricher(gene_L, TERM2GENE = canonical_pathways,minGSSize=min.sz,maxGSSize=max.sz,
                            pvalueCutoff = 1,qvalueCutoff = 1)
          if (is.null(enrich_L)) {
            result_kendall_ora_I=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_kendall_ora_I=enrich_L@result
            result_kendall_ora_I=result_kendall_ora_I[,c("Description","pvalue")]
          }
          colnames(result_kendall_ora_I)=c("geneset_id","P.Value_L")
          result_kendall_ora_I$geneset_id=gsub("\\.", "'", result_kendall_ora_I$geneset_id)
          colnames(result_kendall_ora_A)[2]=paste0(vector[i],".kendall_ora")
          colnames(result_kendall_ora_I)[2]=paste0(vector[i],".kendall_ora")
        }
      }
      if (any(c("spearman_fgsea","spearman_ora") %in% enrichment_method)){
        exp_fGSEA=spearman_test(geneset=exp,value=value)
        if("spearman_fgsea" %in% enrichment_method){
          print("spearman_fgsea will start.")
          alldiff=exp_fGSEA[,c("geneset_id","P.Value","r")]
          min_non_zero=min(alldiff$P.Value[alldiff$P.Value > 0], na.rm = TRUE)
          alldiff$P.Value[alldiff$P.Value == 0]=min_non_zero * 0.1
          alldiff$id=-log10(alldiff$P.Value)*alldiff$r
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
          result_spearman_fgsea=result[,c("pathway","pval","ES")]
          colnames(result_spearman_fgsea)=c("geneset_id","P.Value","ES")
          result_spearman_fgsea$geneset_id=gsub("\\.", "'", result_spearman_fgsea$geneset_id)
          result_spearman_fgsea=na.omit(result_spearman_fgsea)
          result_spearman_fgsea$two=TRUE
          result_spearman_fgsea$invert_A=ifelse(result_spearman_fgsea$ES<0,TRUE,FALSE)
          result_spearman_fgsea$invert_I=ifelse(result_spearman_fgsea$ES>0,TRUE,FALSE)
          result_spearman_fgsea$P.Value.A=two2one(result_spearman_fgsea$P.Value, two = as.logical(result_spearman_fgsea$two), invert = as.logical(result_spearman_fgsea$invert_A))
          result_spearman_fgsea$P.Value.I=two2one(result_spearman_fgsea$P.Value, two = as.logical(result_spearman_fgsea$two), invert = as.logical(result_spearman_fgsea$invert_I))
          result_spearman_fgsea_A=result_spearman_fgsea[,c("geneset_id","P.Value.A")]
          result_spearman_fgsea_I=result_spearman_fgsea[,c("geneset_id","P.Value.I")]
          colnames(result_spearman_fgsea_A)[2]=paste0(vector[i],".spearman_fgsea")
          colnames(result_spearman_fgsea_I)[2]=paste0(vector[i],".spearman_fgsea")
          rm(result_spearman_fgsea)
        }
        if("spearman_ora" %in% enrichment_method){
          print("spearman_ora will start.")
          result_UP=exp_fGSEA[exp_fGSEA$r>=0,]
          result_DN=exp_fGSEA[exp_fGSEA$r<0,]
          result_UP=result_UP[order(result_UP$P.Value, -result_UP$r), ]
          result_DN=result_DN[order(result_DN$P.Value, result_DN$r), ]
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
            result_spearman_ora_A=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_spearman_ora_A=enrich_H@result
            result_spearman_ora_A=result_spearman_ora_A[,c("Description","pvalue")]
          }
          colnames(result_spearman_ora_A)=c("geneset_id","P.Value_H")
          result_spearman_ora_A$geneset_id=gsub("\\.", "'", result_spearman_ora_A$geneset_id)
          enrich_L=enricher(gene_L, TERM2GENE = canonical_pathways,minGSSize=min.sz,maxGSSize=max.sz,
                            pvalueCutoff = 1,qvalueCutoff = 1)
          if (is.null(enrich_L)) {
            result_spearman_ora_I=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_spearman_ora_I=enrich_L@result
            result_spearman_ora_I=result_spearman_ora_I[,c("Description","pvalue")]
          }
          colnames(result_spearman_ora_I)=c("geneset_id","P.Value_L")
          result_spearman_ora_I$geneset_id=gsub("\\.", "'", result_spearman_ora_I$geneset_id)
          colnames(result_spearman_ora_A)[2]=paste0(vector[i],".spearman_ora")
          colnames(result_spearman_ora_I)[2]=paste0(vector[i],".spearman_ora")
        }
      }
      if (any(c("lm_fgsea","lm_ora") %in% enrichment_method)){
        exp_fGSEA=lm_test(geneset=exp,value=value)
        if("lm_fgsea"%in% enrichment_method){
          print("lm_fgsea will start.")
          alldiff=exp_fGSEA[,c("geneset_id","P.Value","r")]
          min_non_zero=min(alldiff$P.Value[alldiff$P.Value > 0], na.rm = TRUE)
          alldiff$P.Value[alldiff$P.Value == 0]=min_non_zero * 0.1
          alldiff$id=-log10(alldiff$P.Value)*alldiff$r
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
          result_lm_fgsea=result[,c("pathway","pval","ES")]
          colnames(result_lm_fgsea)=c("geneset_id","P.Value","ES")
          result_lm_fgsea$geneset_id=gsub("\\.", "'", result_lm_fgsea$geneset_id)
          result_lm_fgsea=na.omit(result_lm_fgsea)
          result_lm_fgsea$two=TRUE
          result_lm_fgsea$invert_A=ifelse(result_lm_fgsea$ES<0,TRUE,FALSE)
          result_lm_fgsea$invert_I=ifelse(result_lm_fgsea$ES>0,TRUE,FALSE)
          result_lm_fgsea$P.Value.A=two2one(result_lm_fgsea$P.Value, two = as.logical(result_lm_fgsea$two), invert = as.logical(result_lm_fgsea$invert_A))
          result_lm_fgsea$P.Value.I=two2one(result_lm_fgsea$P.Value, two = as.logical(result_lm_fgsea$two), invert = as.logical(result_lm_fgsea$invert_I))
          result_lm_fgsea_A=result_lm_fgsea[,c("geneset_id","P.Value.A")]
          result_lm_fgsea_I=result_lm_fgsea[,c("geneset_id","P.Value.I")]
          colnames(result_lm_fgsea_A)[2]=paste0(vector[i],".lm_fgsea")
          colnames(result_lm_fgsea_I)[2]=paste0(vector[i],".lm_fgsea")
          rm(result_lm_fgsea)
        }
        if("lm_ora" %in% enrichment_method){
          print("lm_ora will start.")
          result_UP=exp_fGSEA[exp_fGSEA$r>=0,]
          result_DN=exp_fGSEA[exp_fGSEA$r<0,]
          result_UP=result_UP[order(result_UP$P.Value, -result_UP$r), ]
          result_DN=result_DN[order(result_DN$P.Value, result_DN$r), ]
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
            result_lm_ora_A=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_lm_ora_A=enrich_H@result
            result_lm_ora_A=result_lm_ora_A[,c("Description","pvalue")]
          }
          colnames(result_lm_ora_A)=c("geneset_id","P.Value_H")
          result_lm_ora_A$geneset_id=gsub("\\.", "'", result_lm_ora_A$geneset_id)
          enrich_L=enricher(gene_L, TERM2GENE = canonical_pathways,minGSSize=min.sz,maxGSSize=max.sz,
                            pvalueCutoff = 1,qvalueCutoff = 1)
          if (is.null(enrich_L)) {
            result_lm_ora_I=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_lm_ora_I=enrich_L@result
            result_lm_ora_I=result_lm_ora_I[,c("Description","pvalue")]
          }
          colnames(result_lm_ora_I)=c("geneset_id","P.Value_L")
          result_lm_ora_I$geneset_id=gsub("\\.", "'", result_lm_ora_I$geneset_id)
          colnames(result_lm_ora_A)[2]=paste0(vector[i],".lm_ora")
          colnames(result_lm_ora_I)[2]=paste0(vector[i],".lm_ora")
        }
      }
      if (any(c("biweight_fgsea","biweight_ora") %in% enrichment_method)){
        exp_fGSEA=biweight_test(geneset=exp,value=value)
        if("biweight_fgsea" %in% enrichment_method){
          print("biweight_fgsea will start.")
          alldiff=exp_fGSEA[,c("geneset_id","P.Value","r")]
          min_non_zero=min(alldiff$P.Value[alldiff$P.Value > 0], na.rm = TRUE)
          alldiff$P.Value[alldiff$P.Value == 0]=min_non_zero * 0.1
          alldiff$id=-log10(alldiff$P.Value)*alldiff$r
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
          result_biweight_fgsea=result[,c("pathway","pval","ES")]
          colnames(result_biweight_fgsea)=c("geneset_id","P.Value","ES")
          result_biweight_fgsea$geneset_id=gsub("\\.", "'", result_biweight_fgsea$geneset_id)
          result_biweight_fgsea=na.omit(result_biweight_fgsea)
          result_biweight_fgsea$two=TRUE
          result_biweight_fgsea$invert_A=ifelse(result_biweight_fgsea$ES<0,TRUE,FALSE)
          result_biweight_fgsea$invert_I=ifelse(result_biweight_fgsea$ES>0,TRUE,FALSE)
          result_biweight_fgsea$P.Value.A=two2one(result_biweight_fgsea$P.Value, two = as.logical(result_biweight_fgsea$two), invert = as.logical(result_biweight_fgsea$invert_A))
          result_biweight_fgsea$P.Value.I=two2one(result_biweight_fgsea$P.Value, two = as.logical(result_biweight_fgsea$two), invert = as.logical(result_biweight_fgsea$invert_I))
          result_biweight_fgsea_A=result_biweight_fgsea[,c("geneset_id","P.Value.A")]
          result_biweight_fgsea_I=result_biweight_fgsea[,c("geneset_id","P.Value.I")]
          colnames(result_biweight_fgsea_A)[2]=paste0(vector[i],".biweight_fgsea")
          colnames(result_biweight_fgsea_I)[2]=paste0(vector[i],".biweight_fgsea")
          rm(result_biweight_fgsea)
        }
        if("biweight_ora" %in% enrichment_method){
          print("biweight_ora will start.")
          result_UP=exp_fGSEA[exp_fGSEA$r>=0,]
          result_DN=exp_fGSEA[exp_fGSEA$r<0,]
          result_UP=result_UP[order(result_UP$P.Value, -result_UP$r), ]
          result_DN=result_DN[order(result_DN$P.Value, result_DN$r), ]
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
            result_biweight_ora_A=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_biweight_ora_A=enrich_H@result
            result_biweight_ora_A=result_biweight_ora_A[,c("Description","pvalue")]
          }
          colnames(result_biweight_ora_A)=c("geneset_id","P.Value_H")
          result_biweight_ora_A$geneset_id=gsub("\\.", "'", result_biweight_ora_A$geneset_id)
          enrich_L=enricher(gene_L, TERM2GENE = canonical_pathways,minGSSize=min.sz,maxGSSize=max.sz,
                            pvalueCutoff = 1,qvalueCutoff = 1)
          if (is.null(enrich_L)) {
            result_biweight_ora_I=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_biweight_ora_I=enrich_L@result
            result_biweight_ora_I=result_biweight_ora_I[,c("Description","pvalue")]
          }
          colnames(result_biweight_ora_I)=c("geneset_id","P.Value_L")
          result_biweight_ora_I$geneset_id=gsub("\\.", "'", result_biweight_ora_I$geneset_id)
          colnames(result_biweight_ora_A)[2]=paste0(vector[i],".biweight_ora")
          colnames(result_biweight_ora_I)[2]=paste0(vector[i],".biweight_ora")
        }
      }
      if (any(c("distance_fgsea","distance_ora") %in% enrichment_method)){
        exp_fGSEA=distance_test(geneset=exp,value=value)
        if("distance_fgsea" %in% enrichment_method){
          print("distance_fgsea will start.")
          alldiff=exp_fGSEA[,c("geneset_id","P.Value","r")]
          min_non_zero=min(alldiff$P.Value[alldiff$P.Value > 0], na.rm = TRUE)
          alldiff$P.Value[alldiff$P.Value == 0]=min_non_zero * 0.1
          alldiff$id=-log10(alldiff$P.Value)*alldiff$r
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
          result_distance_fgsea=result[,c("pathway","pval","ES")]
          colnames(result_distance_fgsea)=c("geneset_id","P.Value","ES")
          result_distance_fgsea$geneset_id=gsub("\\.", "'", result_distance_fgsea$geneset_id)
          result_distance_fgsea=na.omit(result_distance_fgsea)
          result_distance_fgsea$two=TRUE
          result_distance_fgsea$invert_A=ifelse(result_distance_fgsea$ES<0,TRUE,FALSE)
          result_distance_fgsea$invert_I=ifelse(result_distance_fgsea$ES>0,TRUE,FALSE)
          result_distance_fgsea$P.Value.A=two2one(result_distance_fgsea$P.Value, two = as.logical(result_distance_fgsea$two), invert = as.logical(result_distance_fgsea$invert_A))
          result_distance_fgsea$P.Value.I=two2one(result_distance_fgsea$P.Value, two = as.logical(result_distance_fgsea$two), invert = as.logical(result_distance_fgsea$invert_I))
          result_distance_fgsea_A=result_distance_fgsea[,c("geneset_id","P.Value.A")]
          result_distance_fgsea_I=result_distance_fgsea[,c("geneset_id","P.Value.I")]
          colnames(result_distance_fgsea_A)[2]=paste0(vector[i],".distance_fgsea")
          colnames(result_distance_fgsea_I)[2]=paste0(vector[i],".distance_fgsea")
          rm(result_distance_fgsea)
        }
        if("distance_ora" %in% enrichment_method){
          print("distance_ora will start.")
          result_UP=exp_fGSEA[exp_fGSEA$r>=0,]
          result_DN=exp_fGSEA[exp_fGSEA$r<0,]
          result_UP=result_UP[order(result_UP$P.Value, -result_UP$r), ]
          result_DN=result_DN[order(result_DN$P.Value, result_DN$r), ]
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
            result_distance_ora_A=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_distance_ora_A=enrich_H@result
            result_distance_ora_A=result_distance_ora_A[,c("Description","pvalue")]
          }
          colnames(result_distance_ora_A)=c("geneset_id","P.Value_H")
          result_distance_ora_A$geneset_id=gsub("\\.", "'", result_distance_ora_A$geneset_id)
          enrich_L=enricher(gene_L, TERM2GENE = canonical_pathways,minGSSize=min.sz,maxGSSize=max.sz,
                            pvalueCutoff = 1,qvalueCutoff = 1)
          if (is.null(enrich_L)) {
            result_distance_ora_I=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_distance_ora_I=enrich_L@result
            result_distance_ora_I=result_distance_ora_I[,c("Description","pvalue")]
          }
          colnames(result_distance_ora_I)=c("geneset_id","P.Value_L")
          result_distance_ora_I$geneset_id=gsub("\\.", "'", result_distance_ora_I$geneset_id)
          colnames(result_distance_ora_A)[2]=paste0(vector[i],".distance_ora")
          colnames(result_distance_ora_I)[2]=paste0(vector[i],".distance_ora")
        }
      }
      if (any(c("percentage_fgsea","percentage_ora") %in% enrichment_method)){
        exp_fGSEA=percentage_test(geneset=exp,value=value)
        if("percentage_fgsea" %in% enrichment_method){
          print("percentage_fgsea will start.")
          alldiff=exp_fGSEA[,c("geneset_id","P.Value","r")]
          min_non_zero=min(alldiff$P.Value[alldiff$P.Value > 0], na.rm = TRUE)
          alldiff$P.Value[alldiff$P.Value == 0]=min_non_zero * 0.1
          alldiff$id=-log10(alldiff$P.Value)*alldiff$r
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
          result_percentage_fgsea=result[,c("pathway","pval","ES")]
          colnames(result_percentage_fgsea)=c("geneset_id","P.Value","ES")
          result_percentage_fgsea$geneset_id=gsub("\\.", "'", result_percentage_fgsea$geneset_id)
          result_percentage_fgsea=na.omit(result_percentage_fgsea)
          result_percentage_fgsea$two=TRUE
          result_percentage_fgsea$invert_A=ifelse(result_percentage_fgsea$ES<0,TRUE,FALSE)
          result_percentage_fgsea$invert_I=ifelse(result_percentage_fgsea$ES>0,TRUE,FALSE)
          result_percentage_fgsea$P.Value.A=two2one(result_percentage_fgsea$P.Value, two = as.logical(result_percentage_fgsea$two), invert = as.logical(result_percentage_fgsea$invert_A))
          result_percentage_fgsea$P.Value.I=two2one(result_percentage_fgsea$P.Value, two = as.logical(result_percentage_fgsea$two), invert = as.logical(result_percentage_fgsea$invert_I))
          result_percentage_fgsea_A=result_percentage_fgsea[,c("geneset_id","P.Value.A")]
          result_percentage_fgsea_I=result_percentage_fgsea[,c("geneset_id","P.Value.I")]
          colnames(result_percentage_fgsea_A)[2]=paste0(vector[i],".percentage_fgsea")
          colnames(result_percentage_fgsea_I)[2]=paste0(vector[i],".percentage_fgsea")
          rm(result_percentage_fgsea)
        }
        if("percentage_ora" %in% enrichment_method){
          print("percentage_ora will start.")
          result_UP=exp_fGSEA[exp_fGSEA$r>=0,]
          result_DN=exp_fGSEA[exp_fGSEA$r<0,]
          result_UP=result_UP[order(result_UP$P.Value, -result_UP$r), ]
          result_DN=result_DN[order(result_DN$P.Value, result_DN$r), ]
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
            result_percentage_ora_A=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_percentage_ora_A=enrich_H@result
            result_percentage_ora_A=result_percentage_ora_A[,c("Description","pvalue")]
          }
          colnames(result_percentage_ora_A)=c("geneset_id","P.Value_H")
          result_percentage_ora_A$geneset_id=gsub("\\.", "'", result_percentage_ora_A$geneset_id)
          enrich_L=enricher(gene_L, TERM2GENE = canonical_pathways,minGSSize=min.sz,maxGSSize=max.sz,
                            pvalueCutoff = 1,qvalueCutoff = 1)
          if (is.null(enrich_L)) {
            result_percentage_ora_I=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_percentage_ora_I=enrich_L@result
            result_percentage_ora_I=result_percentage_ora_I[,c("Description","pvalue")]
          }
          colnames(result_percentage_ora_I)=c("geneset_id","P.Value_L")
          result_percentage_ora_I$geneset_id=gsub("\\.", "'", result_percentage_ora_I$geneset_id)
          colnames(result_percentage_ora_A)[2]=paste0(vector[i],".percentage_ora")
          colnames(result_percentage_ora_I)[2]=paste0(vector[i],".percentage_ora")
        }
      }
      if (any(c("blomqvist_fgsea","blomqvist_ora") %in% enrichment_method)){
        exp_fGSEA=blomqvist_test(geneset=exp,value=value)
        if("blomqvist_fgsea" %in% enrichment_method){
          print("blomqvist_fgsea will start.")
          alldiff=exp_fGSEA[,c("geneset_id","P.Value","r")]
          min_non_zero=min(alldiff$P.Value[alldiff$P.Value > 0], na.rm = TRUE)
          alldiff$P.Value[alldiff$P.Value == 0]=min_non_zero * 0.1
          alldiff$id=-log10(alldiff$P.Value)*alldiff$r
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
          result_blomqvist_fgsea=result[,c("pathway","pval","ES")]
          colnames(result_blomqvist_fgsea)=c("geneset_id","P.Value","ES")
          result_blomqvist_fgsea$geneset_id=gsub("\\.", "'", result_blomqvist_fgsea$geneset_id)
          result_blomqvist_fgsea=na.omit(result_blomqvist_fgsea)
          result_blomqvist_fgsea$two=TRUE
          result_blomqvist_fgsea$invert_A=ifelse(result_blomqvist_fgsea$ES<0,TRUE,FALSE)
          result_blomqvist_fgsea$invert_I=ifelse(result_blomqvist_fgsea$ES>0,TRUE,FALSE)
          result_blomqvist_fgsea$P.Value.A=two2one(result_blomqvist_fgsea$P.Value, two = as.logical(result_blomqvist_fgsea$two), invert = as.logical(result_blomqvist_fgsea$invert_A))
          result_blomqvist_fgsea$P.Value.I=two2one(result_blomqvist_fgsea$P.Value, two = as.logical(result_blomqvist_fgsea$two), invert = as.logical(result_blomqvist_fgsea$invert_I))
          result_blomqvist_fgsea_A=result_blomqvist_fgsea[,c("geneset_id","P.Value.A")]
          result_blomqvist_fgsea_I=result_blomqvist_fgsea[,c("geneset_id","P.Value.I")]
          colnames(result_blomqvist_fgsea_A)[2]=paste0(vector[i],".blomqvist_fgsea")
          colnames(result_blomqvist_fgsea_I)[2]=paste0(vector[i],".blomqvist_fgsea")
          rm(result_blomqvist_fgsea)
        }
        if("blomqvist_ora" %in% enrichment_method){
          print("blomqvist_ora will start.")
          result_UP=exp_fGSEA[exp_fGSEA$r>=0,]
          result_DN=exp_fGSEA[exp_fGSEA$r<0,]
          result_UP=result_UP[order(result_UP$P.Value, -result_UP$r), ]
          result_DN=result_DN[order(result_DN$P.Value, result_DN$r), ]
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
            result_blomqvist_ora_A=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_blomqvist_ora_A=enrich_H@result
            result_blomqvist_ora_A=result_blomqvist_ora_A[,c("Description","pvalue")]
          }
          colnames(result_blomqvist_ora_A)=c("geneset_id","P.Value_H")
          result_blomqvist_ora_A$geneset_id=gsub("\\.", "'", result_blomqvist_ora_A$geneset_id)
          enrich_L=enricher(gene_L, TERM2GENE = canonical_pathways,minGSSize=min.sz,maxGSSize=max.sz,
                            pvalueCutoff = 1,qvalueCutoff = 1)
          if (is.null(enrich_L)) {
            result_blomqvist_ora_I=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_blomqvist_ora_I=enrich_L@result
            result_blomqvist_ora_I=result_blomqvist_ora_I[,c("Description","pvalue")]
          }
          colnames(result_blomqvist_ora_I)=c("geneset_id","P.Value_L")
          result_blomqvist_ora_I$geneset_id=gsub("\\.", "'", result_blomqvist_ora_I$geneset_id)
          colnames(result_blomqvist_ora_A)[2]=paste0(vector[i],".blomqvist_ora")
          colnames(result_blomqvist_ora_I)[2]=paste0(vector[i],".blomqvist_ora")
        }
      }
      if (any(c("hoeffding_fgsea","hoeffding_ora") %in% enrichment_method)){
        exp_fGSEA=hoeffding_test(geneset=exp,value=value)
        if("hoeffding_fgsea" %in% enrichment_method){
          print("hoeffding_fgsea will start.")
          alldiff=exp_fGSEA[,c("geneset_id","P.Value","r")]
          min_non_zero=min(alldiff$P.Value[alldiff$P.Value > 0], na.rm = TRUE)
          alldiff$P.Value[alldiff$P.Value == 0]=min_non_zero * 0.1
          alldiff$id=-log10(alldiff$P.Value)*alldiff$r
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
          result_hoeffding_fgsea=result[,c("pathway","pval","ES")]
          colnames(result_hoeffding_fgsea)=c("geneset_id","P.Value","ES")
          result_hoeffding_fgsea$geneset_id=gsub("\\.", "'", result_hoeffding_fgsea$geneset_id)
          result_hoeffding_fgsea=na.omit(result_hoeffding_fgsea)
          result_hoeffding_fgsea$two=TRUE
          result_hoeffding_fgsea$invert_A=ifelse(result_hoeffding_fgsea$ES<0,TRUE,FALSE)
          result_hoeffding_fgsea$invert_I=ifelse(result_hoeffding_fgsea$ES>0,TRUE,FALSE)
          result_hoeffding_fgsea$P.Value.A=two2one(result_hoeffding_fgsea$P.Value, two = as.logical(result_hoeffding_fgsea$two), invert = as.logical(result_hoeffding_fgsea$invert_A))
          result_hoeffding_fgsea$P.Value.I=two2one(result_hoeffding_fgsea$P.Value, two = as.logical(result_hoeffding_fgsea$two), invert = as.logical(result_hoeffding_fgsea$invert_I))
          result_hoeffding_fgsea_A=result_hoeffding_fgsea[,c("geneset_id","P.Value.A")]
          result_hoeffding_fgsea_I=result_hoeffding_fgsea[,c("geneset_id","P.Value.I")]
          colnames(result_hoeffding_fgsea_A)[2]=paste0(vector[i],".hoeffding_fgsea")
          colnames(result_hoeffding_fgsea_I)[2]=paste0(vector[i],".hoeffding_fgsea")
          rm(result_hoeffding_fgsea)
        }
        if("hoeffding_ora" %in% enrichment_method){
          print("hoeffding_ora will start.")
          result_UP=exp_fGSEA[exp_fGSEA$r>=0,]
          result_DN=exp_fGSEA[exp_fGSEA$r<0,]
          result_UP=result_UP[order(result_UP$P.Value, -result_UP$r), ]
          result_DN=result_DN[order(result_DN$P.Value, result_DN$r), ]
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
            result_hoeffding_ora_A=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_hoeffding_ora_A=enrich_H@result
            result_hoeffding_ora_A=result_hoeffding_ora_A[,c("Description","pvalue")]
          }
          colnames(result_hoeffding_ora_A)=c("geneset_id","P.Value_H")
          result_hoeffding_ora_A$geneset_id=gsub("\\.", "'", result_hoeffding_ora_A$geneset_id)
          enrich_L=enricher(gene_L, TERM2GENE = canonical_pathways,minGSSize=min.sz,maxGSSize=max.sz,
                            pvalueCutoff = 1,qvalueCutoff = 1)
          if (is.null(enrich_L)) {
            result_hoeffding_ora_I=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_hoeffding_ora_I=enrich_L@result
            result_hoeffding_ora_I=result_hoeffding_ora_I[,c("Description","pvalue")]
          }
          colnames(result_hoeffding_ora_I)=c("geneset_id","P.Value_L")
          result_hoeffding_ora_I$geneset_id=gsub("\\.", "'", result_hoeffding_ora_I$geneset_id)
          colnames(result_hoeffding_ora_A)[2]=paste0(vector[i],".hoeffding_ora")
          colnames(result_hoeffding_ora_I)[2]=paste0(vector[i],".hoeffding_ora")
        }
      }
      if (any(c("gamma_fgsea","gamma_ora") %in% enrichment_method)){
        exp_fGSEA=gamma_test(geneset=exp,value=value)
        if("gamma_fgsea" %in% enrichment_method){
          print("gamma_fgsea will start.")
          alldiff=exp_fGSEA[,c("geneset_id","P.Value","r")]
          min_non_zero=min(alldiff$P.Value[alldiff$P.Value > 0], na.rm = TRUE)
          alldiff$P.Value[alldiff$P.Value == 0]=min_non_zero * 0.1
          alldiff$id=-log10(alldiff$P.Value)*alldiff$r
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
          result_gamma_fgsea=result[,c("pathway","pval","ES")]
          colnames(result_gamma_fgsea)=c("geneset_id","P.Value","ES")
          result_gamma_fgsea$geneset_id=gsub("\\.", "'", result_gamma_fgsea$geneset_id)
          result_gamma_fgsea=na.omit(result_gamma_fgsea)
          result_gamma_fgsea$two=TRUE
          result_gamma_fgsea$invert_A=ifelse(result_gamma_fgsea$ES<0,TRUE,FALSE)
          result_gamma_fgsea$invert_I=ifelse(result_gamma_fgsea$ES>0,TRUE,FALSE)
          result_gamma_fgsea$P.Value.A=two2one(result_gamma_fgsea$P.Value, two = as.logical(result_gamma_fgsea$two), invert = as.logical(result_gamma_fgsea$invert_A))
          result_gamma_fgsea$P.Value.I=two2one(result_gamma_fgsea$P.Value, two = as.logical(result_gamma_fgsea$two), invert = as.logical(result_gamma_fgsea$invert_I))
          result_gamma_fgsea_A=result_gamma_fgsea[,c("geneset_id","P.Value.A")]
          result_gamma_fgsea_I=result_gamma_fgsea[,c("geneset_id","P.Value.I")]
          colnames(result_gamma_fgsea_A)[2]=paste0(vector[i],".gamma_fgsea")
          colnames(result_gamma_fgsea_I)[2]=paste0(vector[i],".gamma_fgsea")
          rm(result_gamma_fgsea)
        }
        if("gamma_ora" %in% enrichment_method){
          print("gamma_ora will start.")
          result_UP=exp_fGSEA[exp_fGSEA$r>=0,]
          result_DN=exp_fGSEA[exp_fGSEA$r<0,]
          result_UP=result_UP[order(result_UP$P.Value, -result_UP$r), ]
          result_DN=result_DN[order(result_DN$P.Value, result_DN$r), ]
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
            result_gamma_ora_A=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_gamma_ora_A=enrich_H@result
            result_gamma_ora_A=result_gamma_ora_A[,c("Description","pvalue")]
          }
          colnames(result_gamma_ora_A)=c("geneset_id","P.Value_H")
          result_gamma_ora_A$geneset_id=gsub("\\.", "'", result_gamma_ora_A$geneset_id)
          enrich_L=enricher(gene_L, TERM2GENE = canonical_pathways,minGSSize=min.sz,maxGSSize=max.sz,
                            pvalueCutoff = 1,qvalueCutoff = 1)
          if (is.null(enrich_L)) {
            result_gamma_ora_I=data.frame(Description=label_genesets,pvalue=rep(1, length(label_genesets)))
          }else{
            result_gamma_ora_I=enrich_L@result
            result_gamma_ora_I=result_gamma_ora_I[,c("Description","pvalue")]
          }
          colnames(result_gamma_ora_I)=c("geneset_id","P.Value_L")
          result_gamma_ora_I$geneset_id=gsub("\\.", "'", result_gamma_ora_I$geneset_id)
          colnames(result_gamma_ora_A)[2]=paste0(vector[i],".gamma_ora")
          colnames(result_gamma_ora_I)[2]=paste0(vector[i],".gamma_ora")
        }
      }
    } else {
      print("Not using fgsea and/or ora")
    }
    ##########################################consensus
    methods_to_check=c("consensus_pearson","consensus_kendall","consensus_spearman","consensus_lm","consensus_biweight","consensus_distance","consensus_percentage","consensus_blomqvist","consensus_hoeffding","consensus_gamma")
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
      if("consensus_pearson" %in% enrichment_method){
        print("consensus_pearson will start.")
        result_consensus_pearson=pearson_test(geneset=geneset,value=value)
        result_consensus_pearson$two=TRUE
        result_consensus_pearson$invert_A=ifelse(result_consensus_pearson$r<=0,TRUE,FALSE)
        result_consensus_pearson$invert_I=ifelse(result_consensus_pearson$r>0,TRUE,FALSE)
        result_consensus_pearson$P.Value.A=two2one(result_consensus_pearson$P.Value, two = as.logical(result_consensus_pearson$two), invert = as.logical(result_consensus_pearson$invert_A))
        result_consensus_pearson$P.Value.I=two2one(result_consensus_pearson$P.Value, two = as.logical(result_consensus_pearson$two), invert = as.logical(result_consensus_pearson$invert_I))
        result_consensus_pearson_A=result_consensus_pearson[,c("geneset_id","P.Value.A")]
        result_consensus_pearson_I=result_consensus_pearson[,c("geneset_id","P.Value.I")]
        colnames(result_consensus_pearson_A)[2]=paste0(vector[i],".consensus_pearson")
        colnames(result_consensus_pearson_I)[2]=paste0(vector[i],".consensus_pearson")
        rm(result_consensus_pearson)
      }
      if("consensus_kendall" %in% enrichment_method){
        print("consensus_kendall will start.")
        result_consensus_kendall=kendall_test(geneset=geneset,value=value)
        result_consensus_kendall$two=TRUE
        result_consensus_kendall$invert_A=ifelse(result_consensus_kendall$r<=0,TRUE,FALSE)
        result_consensus_kendall$invert_I=ifelse(result_consensus_kendall$r>0,TRUE,FALSE)
        result_consensus_kendall$P.Value.A=two2one(result_consensus_kendall$P.Value, two = as.logical(result_consensus_kendall$two), invert = as.logical(result_consensus_kendall$invert_A))
        result_consensus_kendall$P.Value.I=two2one(result_consensus_kendall$P.Value, two = as.logical(result_consensus_kendall$two), invert = as.logical(result_consensus_kendall$invert_I))
        result_consensus_kendall_A=result_consensus_kendall[,c("geneset_id","P.Value.A")]
        result_consensus_kendall_I=result_consensus_kendall[,c("geneset_id","P.Value.I")]
        colnames(result_consensus_kendall_A)[2]=paste0(vector[i],".consensus_kendall")
        colnames(result_consensus_kendall_I)[2]=paste0(vector[i],".consensus_kendall")
        rm(result_consensus_kendall)
      }
      if("consensus_spearman" %in% enrichment_method){
        print("consensus_spearman will start.")
        result_consensus_spearman=spearman_test(geneset=geneset,value=value)
        result_consensus_spearman$two=TRUE
        result_consensus_spearman$invert_A=ifelse(result_consensus_spearman$r<=0,TRUE,FALSE)
        result_consensus_spearman$invert_I=ifelse(result_consensus_spearman$r>0,TRUE,FALSE)
        result_consensus_spearman$P.Value.A=two2one(result_consensus_spearman$P.Value, two = as.logical(result_consensus_spearman$two), invert = as.logical(result_consensus_spearman$invert_A))
        result_consensus_spearman$P.Value.I=two2one(result_consensus_spearman$P.Value, two = as.logical(result_consensus_spearman$two), invert = as.logical(result_consensus_spearman$invert_I))
        result_consensus_spearman_A=result_consensus_spearman[,c("geneset_id","P.Value.A")]
        result_consensus_spearman_I=result_consensus_spearman[,c("geneset_id","P.Value.I")]
        colnames(result_consensus_spearman_A)[2]=paste0(vector[i],".consensus_spearman")
        colnames(result_consensus_spearman_I)[2]=paste0(vector[i],".consensus_spearman")
        rm(result_consensus_spearman)
      }
      if("consensus_lm" %in% enrichment_method){
        print("consensus_lm will start.")
        result_consensus_lm=lm_test(geneset=geneset,value=value)
        result_consensus_lm$two=TRUE
        result_consensus_lm$invert_A=ifelse(result_consensus_lm$r<=0,TRUE,FALSE)
        result_consensus_lm$invert_I=ifelse(result_consensus_lm$r>0,TRUE,FALSE)
        result_consensus_lm$P.Value.A=two2one(result_consensus_lm$P.Value, two = as.logical(result_consensus_lm$two), invert = as.logical(result_consensus_lm$invert_A))
        result_consensus_lm$P.Value.I=two2one(result_consensus_lm$P.Value, two = as.logical(result_consensus_lm$two), invert = as.logical(result_consensus_lm$invert_I))
        result_consensus_lm_A=result_consensus_lm[,c("geneset_id","P.Value.A")]
        result_consensus_lm_I=result_consensus_lm[,c("geneset_id","P.Value.I")]
        colnames(result_consensus_lm_A)[2]=paste0(vector[i],".consensus_lm")
        colnames(result_consensus_lm_I)[2]=paste0(vector[i],".consensus_lm")
        rm(result_consensus_lm)
      }
      if("consensus_biweight" %in% enrichment_method){
        print("consensus_biweight will start.")
        result_consensus_biweight=biweight_test(geneset=geneset,value=value)
        result_consensus_biweight$two=TRUE
        result_consensus_biweight$invert_A=ifelse(result_consensus_biweight$r<=0,TRUE,FALSE)
        result_consensus_biweight$invert_I=ifelse(result_consensus_biweight$r>0,TRUE,FALSE)
        result_consensus_biweight$P.Value.A=two2one(result_consensus_biweight$P.Value, two = as.logical(result_consensus_biweight$two), invert = as.logical(result_consensus_biweight$invert_A))
        result_consensus_biweight$P.Value.I=two2one(result_consensus_biweight$P.Value, two = as.logical(result_consensus_biweight$two), invert = as.logical(result_consensus_biweight$invert_I))
        result_consensus_biweight_A=result_consensus_biweight[,c("geneset_id","P.Value.A")]
        result_consensus_biweight_I=result_consensus_biweight[,c("geneset_id","P.Value.I")]
        colnames(result_consensus_biweight_A)[2]=paste0(vector[i],".consensus_biweight")
        colnames(result_consensus_biweight_I)[2]=paste0(vector[i],".consensus_biweight")
        rm(result_consensus_biweight)
      }
      if("consensus_distance" %in% enrichment_method){
        print("consensus_distance will start.")
        result_consensus_distance=distance_test(geneset=geneset,value=value)
        result_consensus_distance$two=TRUE
        result_consensus_distance$invert_A=ifelse(result_consensus_distance$r<=0,TRUE,FALSE)
        result_consensus_distance$invert_I=ifelse(result_consensus_distance$r>0,TRUE,FALSE)
        result_consensus_distance$P.Value.A=two2one(result_consensus_distance$P.Value, two = as.logical(result_consensus_distance$two), invert = as.logical(result_consensus_distance$invert_A))
        result_consensus_distance$P.Value.I=two2one(result_consensus_distance$P.Value, two = as.logical(result_consensus_distance$two), invert = as.logical(result_consensus_distance$invert_I))
        result_consensus_distance_A=result_consensus_distance[,c("geneset_id","P.Value.A")]
        result_consensus_distance_I=result_consensus_distance[,c("geneset_id","P.Value.I")]
        colnames(result_consensus_distance_A)[2]=paste0(vector[i],".consensus_distance")
        colnames(result_consensus_distance_I)[2]=paste0(vector[i],".consensus_distance")
        rm(result_consensus_distance)
      }
      if("consensus_percentage" %in% enrichment_method){
        print("consensus_percentage will start.")
        result_consensus_percentage=percentage_test(geneset=geneset,value=value)
        result_consensus_percentage$two=TRUE
        result_consensus_percentage$invert_A=ifelse(result_consensus_percentage$r<=0,TRUE,FALSE)
        result_consensus_percentage$invert_I=ifelse(result_consensus_percentage$r>0,TRUE,FALSE)
        result_consensus_percentage$P.Value.A=two2one(result_consensus_percentage$P.Value, two = as.logical(result_consensus_percentage$two), invert = as.logical(result_consensus_percentage$invert_A))
        result_consensus_percentage$P.Value.I=two2one(result_consensus_percentage$P.Value, two = as.logical(result_consensus_percentage$two), invert = as.logical(result_consensus_percentage$invert_I))
        result_consensus_percentage_A=result_consensus_percentage[,c("geneset_id","P.Value.A")]
        result_consensus_percentage_I=result_consensus_percentage[,c("geneset_id","P.Value.I")]
        colnames(result_consensus_percentage_A)[2]=paste0(vector[i],".consensus_percentage")
        colnames(result_consensus_percentage_I)[2]=paste0(vector[i],".consensus_percentage")
        rm(result_consensus_percentage)
      }
      if("consensus_blomqvist" %in% enrichment_method){
        print("consensus_blomqvist will start.")
        result_consensus_blomqvist=blomqvist_test(geneset=geneset,value=value)
        result_consensus_blomqvist$two=TRUE
        result_consensus_blomqvist$invert_A=ifelse(result_consensus_blomqvist$r<=0,TRUE,FALSE)
        result_consensus_blomqvist$invert_I=ifelse(result_consensus_blomqvist$r>0,TRUE,FALSE)
        result_consensus_blomqvist$P.Value.A=two2one(result_consensus_blomqvist$P.Value, two = as.logical(result_consensus_blomqvist$two), invert = as.logical(result_consensus_blomqvist$invert_A))
        result_consensus_blomqvist$P.Value.I=two2one(result_consensus_blomqvist$P.Value, two = as.logical(result_consensus_blomqvist$two), invert = as.logical(result_consensus_blomqvist$invert_I))
        result_consensus_blomqvist_A=result_consensus_blomqvist[,c("geneset_id","P.Value.A")]
        result_consensus_blomqvist_I=result_consensus_blomqvist[,c("geneset_id","P.Value.I")]
        colnames(result_consensus_blomqvist_A)[2]=paste0(vector[i],".consensus_blomqvist")
        colnames(result_consensus_blomqvist_I)[2]=paste0(vector[i],".consensus_blomqvist")
        rm(result_consensus_blomqvist)
      }
      if("consensus_hoeffding" %in% enrichment_method){
        print("consensus_hoeffding will start.")
        result_consensus_hoeffding=hoeffding_test(geneset=geneset,value=value)
        result_consensus_hoeffding$two=TRUE
        result_consensus_hoeffding$invert_A=ifelse(result_consensus_hoeffding$r<=0,TRUE,FALSE)
        result_consensus_hoeffding$invert_I=ifelse(result_consensus_hoeffding$r>0,TRUE,FALSE)
        result_consensus_hoeffding$P.Value.A=two2one(result_consensus_hoeffding$P.Value, two = as.logical(result_consensus_hoeffding$two), invert = as.logical(result_consensus_hoeffding$invert_A))
        result_consensus_hoeffding$P.Value.I=two2one(result_consensus_hoeffding$P.Value, two = as.logical(result_consensus_hoeffding$two), invert = as.logical(result_consensus_hoeffding$invert_I))
        result_consensus_hoeffding_A=result_consensus_hoeffding[,c("geneset_id","P.Value.A")]
        result_consensus_hoeffding_I=result_consensus_hoeffding[,c("geneset_id","P.Value.I")]
        colnames(result_consensus_hoeffding_A)[2]=paste0(vector[i],".consensus_hoeffding")
        colnames(result_consensus_hoeffding_I)[2]=paste0(vector[i],".consensus_hoeffding")
        rm(result_consensus_hoeffding)
      }
      if("consensus_gamma" %in% enrichment_method){
        print("consensus_gamma will start.")
        result_consensus_gamma=gamma_test(geneset=geneset,value=value)
        result_consensus_gamma$two=TRUE
        result_consensus_gamma$invert_A=ifelse(result_consensus_gamma$r<=0,TRUE,FALSE)
        result_consensus_gamma$invert_I=ifelse(result_consensus_gamma$r>0,TRUE,FALSE)
        result_consensus_gamma$P.Value.A=two2one(result_consensus_gamma$P.Value, two = as.logical(result_consensus_gamma$two), invert = as.logical(result_consensus_gamma$invert_A))
        result_consensus_gamma$P.Value.I=two2one(result_consensus_gamma$P.Value, two = as.logical(result_consensus_gamma$two), invert = as.logical(result_consensus_gamma$invert_I))
        result_consensus_gamma_A=result_consensus_gamma[,c("geneset_id","P.Value.A")]
        result_consensus_gamma_I=result_consensus_gamma[,c("geneset_id","P.Value.I")]
        colnames(result_consensus_gamma_A)[2]=paste0(vector[i],".consensus_gamma")
        colnames(result_consensus_gamma_I)[2]=paste0(vector[i],".consensus_gamma")
        rm(result_consensus_gamma)
      }
    }else {
      print("Not using consensus.")
    }
    ##########################################mdt
    methods_to_check=c("mdt_pearson","mdt_kendall","mdt_spearman","mdt_lm","mdt_biweight","mdt_distance","mdt_percentage","mdt_blomqvist","mdt_hoeffding","mdt_gamma")
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
      if("mdt_pearson" %in% enrichment_method){
        print("mdt_pearson will start.")
        result_mdt_pearson=pearson_test(geneset=geneset,value=value)
        result_mdt_pearson$two=TRUE
        result_mdt_pearson$invert_A=ifelse(result_mdt_pearson$r<=0,TRUE,FALSE)
        result_mdt_pearson$invert_I=ifelse(result_mdt_pearson$r>0,TRUE,FALSE)
        result_mdt_pearson$P.Value.A=two2one(result_mdt_pearson$P.Value, two = as.logical(result_mdt_pearson$two), invert = as.logical(result_mdt_pearson$invert_A))
        result_mdt_pearson$P.Value.I=two2one(result_mdt_pearson$P.Value, two = as.logical(result_mdt_pearson$two), invert = as.logical(result_mdt_pearson$invert_I))
        result_mdt_pearson_A=result_mdt_pearson[,c("geneset_id","P.Value.A")]
        result_mdt_pearson_I=result_mdt_pearson[,c("geneset_id","P.Value.I")]
        colnames(result_mdt_pearson_A)[2]=paste0(vector[i],".mdt_pearson")
        colnames(result_mdt_pearson_I)[2]=paste0(vector[i],".mdt_pearson")
        rm(result_mdt_pearson)
      }
      if("mdt_kendall" %in% enrichment_method){
        print("mdt_kendall will start.")
        result_mdt_kendall=kendall_test(geneset=geneset,value=value)
        result_mdt_kendall$two=TRUE
        result_mdt_kendall$invert_A=ifelse(result_mdt_kendall$r<=0,TRUE,FALSE)
        result_mdt_kendall$invert_I=ifelse(result_mdt_kendall$r>0,TRUE,FALSE)
        result_mdt_kendall$P.Value.A=two2one(result_mdt_kendall$P.Value, two = as.logical(result_mdt_kendall$two), invert = as.logical(result_mdt_kendall$invert_A))
        result_mdt_kendall$P.Value.I=two2one(result_mdt_kendall$P.Value, two = as.logical(result_mdt_kendall$two), invert = as.logical(result_mdt_kendall$invert_I))
        result_mdt_kendall_A=result_mdt_kendall[,c("geneset_id","P.Value.A")]
        result_mdt_kendall_I=result_mdt_kendall[,c("geneset_id","P.Value.I")]
        colnames(result_mdt_kendall_A)[2]=paste0(vector[i],".mdt_kendall")
        colnames(result_mdt_kendall_I)[2]=paste0(vector[i],".mdt_kendall")
        rm(result_mdt_kendall)
      }
      if("mdt_spearman" %in% enrichment_method){
        print("mdt_spearman will start.")
        result_mdt_spearman=spearman_test(geneset=geneset,value=value)
        result_mdt_spearman$two=TRUE
        result_mdt_spearman$invert_A=ifelse(result_mdt_spearman$r<=0,TRUE,FALSE)
        result_mdt_spearman$invert_I=ifelse(result_mdt_spearman$r>0,TRUE,FALSE)
        result_mdt_spearman$P.Value.A=two2one(result_mdt_spearman$P.Value, two = as.logical(result_mdt_spearman$two), invert = as.logical(result_mdt_spearman$invert_A))
        result_mdt_spearman$P.Value.I=two2one(result_mdt_spearman$P.Value, two = as.logical(result_mdt_spearman$two), invert = as.logical(result_mdt_spearman$invert_I))
        result_mdt_spearman_A=result_mdt_spearman[,c("geneset_id","P.Value.A")]
        result_mdt_spearman_I=result_mdt_spearman[,c("geneset_id","P.Value.I")]
        colnames(result_mdt_spearman_A)[2]=paste0(vector[i],".mdt_spearman")
        colnames(result_mdt_spearman_I)[2]=paste0(vector[i],".mdt_spearman")
        rm(result_mdt_spearman)
      }
      if("mdt_lm" %in% enrichment_method){
        print("mdt_lm will start.")
        result_mdt_lm=lm_test(geneset=geneset,value=value)
        result_mdt_lm$two=TRUE
        result_mdt_lm$invert_A=ifelse(result_mdt_lm$r<=0,TRUE,FALSE)
        result_mdt_lm$invert_I=ifelse(result_mdt_lm$r>0,TRUE,FALSE)
        result_mdt_lm$P.Value.A=two2one(result_mdt_lm$P.Value, two = as.logical(result_mdt_lm$two), invert = as.logical(result_mdt_lm$invert_A))
        result_mdt_lm$P.Value.I=two2one(result_mdt_lm$P.Value, two = as.logical(result_mdt_lm$two), invert = as.logical(result_mdt_lm$invert_I))
        result_mdt_lm_A=result_mdt_lm[,c("geneset_id","P.Value.A")]
        result_mdt_lm_I=result_mdt_lm[,c("geneset_id","P.Value.I")]
        colnames(result_mdt_lm_A)[2]=paste0(vector[i],".mdt_lm")
        colnames(result_mdt_lm_I)[2]=paste0(vector[i],".mdt_lm")
        rm(result_mdt_lm)
      }
      if("mdt_biweight" %in% enrichment_method){
        print("mdt_biweight will start.")
        result_mdt_biweight=biweight_test(geneset=geneset,value=value)
        result_mdt_biweight$two=TRUE
        result_mdt_biweight$invert_A=ifelse(result_mdt_biweight$r<=0,TRUE,FALSE)
        result_mdt_biweight$invert_I=ifelse(result_mdt_biweight$r>0,TRUE,FALSE)
        result_mdt_biweight$P.Value.A=two2one(result_mdt_biweight$P.Value, two = as.logical(result_mdt_biweight$two), invert = as.logical(result_mdt_biweight$invert_A))
        result_mdt_biweight$P.Value.I=two2one(result_mdt_biweight$P.Value, two = as.logical(result_mdt_biweight$two), invert = as.logical(result_mdt_biweight$invert_I))
        result_mdt_biweight_A=result_mdt_biweight[,c("geneset_id","P.Value.A")]
        result_mdt_biweight_I=result_mdt_biweight[,c("geneset_id","P.Value.I")]
        colnames(result_mdt_biweight_A)[2]=paste0(vector[i],".mdt_biweight")
        colnames(result_mdt_biweight_I)[2]=paste0(vector[i],".mdt_biweight")
        rm(result_mdt_biweight)
      }
      if("mdt_distance" %in% enrichment_method){
        print("mdt_distance will start.")
        result_mdt_distance=distance_test(geneset=geneset,value=value)
        result_mdt_distance$two=TRUE
        result_mdt_distance$invert_A=ifelse(result_mdt_distance$r<=0,TRUE,FALSE)
        result_mdt_distance$invert_I=ifelse(result_mdt_distance$r>0,TRUE,FALSE)
        result_mdt_distance$P.Value.A=two2one(result_mdt_distance$P.Value, two = as.logical(result_mdt_distance$two), invert = as.logical(result_mdt_distance$invert_A))
        result_mdt_distance$P.Value.I=two2one(result_mdt_distance$P.Value, two = as.logical(result_mdt_distance$two), invert = as.logical(result_mdt_distance$invert_I))
        result_mdt_distance_A=result_mdt_distance[,c("geneset_id","P.Value.A")]
        result_mdt_distance_I=result_mdt_distance[,c("geneset_id","P.Value.I")]
        colnames(result_mdt_distance_A)[2]=paste0(vector[i],".mdt_distance")
        colnames(result_mdt_distance_I)[2]=paste0(vector[i],".mdt_distance")
        rm(result_mdt_distance)
      }
      if("mdt_percentage" %in% enrichment_method){
        print("mdt_percentage will start.")
        result_mdt_percentage=percentage_test(geneset=geneset,value=value)
        result_mdt_percentage$two=TRUE
        result_mdt_percentage$invert_A=ifelse(result_mdt_percentage$r<=0,TRUE,FALSE)
        result_mdt_percentage$invert_I=ifelse(result_mdt_percentage$r>0,TRUE,FALSE)
        result_mdt_percentage$P.Value.A=two2one(result_mdt_percentage$P.Value, two = as.logical(result_mdt_percentage$two), invert = as.logical(result_mdt_percentage$invert_A))
        result_mdt_percentage$P.Value.I=two2one(result_mdt_percentage$P.Value, two = as.logical(result_mdt_percentage$two), invert = as.logical(result_mdt_percentage$invert_I))
        result_mdt_percentage_A=result_mdt_percentage[,c("geneset_id","P.Value.A")]
        result_mdt_percentage_I=result_mdt_percentage[,c("geneset_id","P.Value.I")]
        colnames(result_mdt_percentage_A)[2]=paste0(vector[i],".mdt_percentage")
        colnames(result_mdt_percentage_I)[2]=paste0(vector[i],".mdt_percentage")
        rm(result_mdt_percentage)
      }
      if("mdt_blomqvist" %in% enrichment_method){
        print("mdt_blomqvist will start.")
        result_mdt_blomqvist=blomqvist_test(geneset=geneset,value=value)
        result_mdt_blomqvist$two=TRUE
        result_mdt_blomqvist$invert_A=ifelse(result_mdt_blomqvist$r<=0,TRUE,FALSE)
        result_mdt_blomqvist$invert_I=ifelse(result_mdt_blomqvist$r>0,TRUE,FALSE)
        result_mdt_blomqvist$P.Value.A=two2one(result_mdt_blomqvist$P.Value, two = as.logical(result_mdt_blomqvist$two), invert = as.logical(result_mdt_blomqvist$invert_A))
        result_mdt_blomqvist$P.Value.I=two2one(result_mdt_blomqvist$P.Value, two = as.logical(result_mdt_blomqvist$two), invert = as.logical(result_mdt_blomqvist$invert_I))
        result_mdt_blomqvist_A=result_mdt_blomqvist[,c("geneset_id","P.Value.A")]
        result_mdt_blomqvist_I=result_mdt_blomqvist[,c("geneset_id","P.Value.I")]
        colnames(result_mdt_blomqvist_A)[2]=paste0(vector[i],".mdt_blomqvist")
        colnames(result_mdt_blomqvist_I)[2]=paste0(vector[i],".mdt_blomqvist")
        rm(result_mdt_blomqvist)
      }
      if("mdt_hoeffding" %in% enrichment_method){
        print("mdt_hoeffding will start.")
        result_mdt_hoeffding=hoeffding_test(geneset=geneset,value=value)
        result_mdt_hoeffding$two=TRUE
        result_mdt_hoeffding$invert_A=ifelse(result_mdt_hoeffding$r<=0,TRUE,FALSE)
        result_mdt_hoeffding$invert_I=ifelse(result_mdt_hoeffding$r>0,TRUE,FALSE)
        result_mdt_hoeffding$P.Value.A=two2one(result_mdt_hoeffding$P.Value, two = as.logical(result_mdt_hoeffding$two), invert = as.logical(result_mdt_hoeffding$invert_A))
        result_mdt_hoeffding$P.Value.I=two2one(result_mdt_hoeffding$P.Value, two = as.logical(result_mdt_hoeffding$two), invert = as.logical(result_mdt_hoeffding$invert_I))
        result_mdt_hoeffding_A=result_mdt_hoeffding[,c("geneset_id","P.Value.A")]
        result_mdt_hoeffding_I=result_mdt_hoeffding[,c("geneset_id","P.Value.I")]
        colnames(result_mdt_hoeffding_A)[2]=paste0(vector[i],".mdt_hoeffding")
        colnames(result_mdt_hoeffding_I)[2]=paste0(vector[i],".mdt_hoeffding")
        rm(result_mdt_hoeffding)
      }
      if("mdt_gamma" %in% enrichment_method){
        print("mdt_gamma will start.")
        result_mdt_gamma=gamma_test(geneset=geneset,value=value)
        result_mdt_gamma$two=TRUE
        result_mdt_gamma$invert_A=ifelse(result_mdt_gamma$r<=0,TRUE,FALSE)
        result_mdt_gamma$invert_I=ifelse(result_mdt_gamma$r>0,TRUE,FALSE)
        result_mdt_gamma$P.Value.A=two2one(result_mdt_gamma$P.Value, two = as.logical(result_mdt_gamma$two), invert = as.logical(result_mdt_gamma$invert_A))
        result_mdt_gamma$P.Value.I=two2one(result_mdt_gamma$P.Value, two = as.logical(result_mdt_gamma$two), invert = as.logical(result_mdt_gamma$invert_I))
        result_mdt_gamma_A=result_mdt_gamma[,c("geneset_id","P.Value.A")]
        result_mdt_gamma_I=result_mdt_gamma[,c("geneset_id","P.Value.I")]
        colnames(result_mdt_gamma_A)[2]=paste0(vector[i],".mdt_gamma")
        colnames(result_mdt_gamma_I)[2]=paste0(vector[i],".mdt_gamma")
        rm(result_mdt_gamma)
      }
    }else {
      print("Not using mdt.")
    }
    ##########################################mlm
    methods_to_check=c("mlm_pearson","mlm_kendall","mlm_spearman","mlm_lm","mlm_biweight","mlm_distance","mlm_percentage","mlm_blomqvist","mlm_hoeffding","mlm_gamma")
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
      if("mlm_pearson" %in% enrichment_method){
        print("mlm_pearson will start.")
        result_mlm_pearson=pearson_test(geneset=geneset,value=value)
        result_mlm_pearson$two=TRUE
        result_mlm_pearson$invert_A=ifelse(result_mlm_pearson$r<=0,TRUE,FALSE)
        result_mlm_pearson$invert_I=ifelse(result_mlm_pearson$r>0,TRUE,FALSE)
        result_mlm_pearson$P.Value.A=two2one(result_mlm_pearson$P.Value, two = as.logical(result_mlm_pearson$two), invert = as.logical(result_mlm_pearson$invert_A))
        result_mlm_pearson$P.Value.I=two2one(result_mlm_pearson$P.Value, two = as.logical(result_mlm_pearson$two), invert = as.logical(result_mlm_pearson$invert_I))
        result_mlm_pearson_A=result_mlm_pearson[,c("geneset_id","P.Value.A")]
        result_mlm_pearson_I=result_mlm_pearson[,c("geneset_id","P.Value.I")]
        colnames(result_mlm_pearson_A)[2]=paste0(vector[i],".mlm_pearson")
        colnames(result_mlm_pearson_I)[2]=paste0(vector[i],".mlm_pearson")
        rm(result_mlm_pearson)
      }
      if("mlm_kendall" %in% enrichment_method){
        print("mlm_kendall will start.")
        result_mlm_kendall=kendall_test(geneset=geneset,value=value)
        result_mlm_kendall$two=TRUE
        result_mlm_kendall$invert_A=ifelse(result_mlm_kendall$r<=0,TRUE,FALSE)
        result_mlm_kendall$invert_I=ifelse(result_mlm_kendall$r>0,TRUE,FALSE)
        result_mlm_kendall$P.Value.A=two2one(result_mlm_kendall$P.Value, two = as.logical(result_mlm_kendall$two), invert = as.logical(result_mlm_kendall$invert_A))
        result_mlm_kendall$P.Value.I=two2one(result_mlm_kendall$P.Value, two = as.logical(result_mlm_kendall$two), invert = as.logical(result_mlm_kendall$invert_I))
        result_mlm_kendall_A=result_mlm_kendall[,c("geneset_id","P.Value.A")]
        result_mlm_kendall_I=result_mlm_kendall[,c("geneset_id","P.Value.I")]
        colnames(result_mlm_kendall_A)[2]=paste0(vector[i],".mlm_kendall")
        colnames(result_mlm_kendall_I)[2]=paste0(vector[i],".mlm_kendall")
        rm(result_mlm_kendall)
      }
      if("mlm_spearman" %in% enrichment_method){
        print("mlm_spearman will start.")
        result_mlm_spearman=spearman_test(geneset=geneset,value=value)
        result_mlm_spearman$two=TRUE
        result_mlm_spearman$invert_A=ifelse(result_mlm_spearman$r<=0,TRUE,FALSE)
        result_mlm_spearman$invert_I=ifelse(result_mlm_spearman$r>0,TRUE,FALSE)
        result_mlm_spearman$P.Value.A=two2one(result_mlm_spearman$P.Value, two = as.logical(result_mlm_spearman$two), invert = as.logical(result_mlm_spearman$invert_A))
        result_mlm_spearman$P.Value.I=two2one(result_mlm_spearman$P.Value, two = as.logical(result_mlm_spearman$two), invert = as.logical(result_mlm_spearman$invert_I))
        result_mlm_spearman_A=result_mlm_spearman[,c("geneset_id","P.Value.A")]
        result_mlm_spearman_I=result_mlm_spearman[,c("geneset_id","P.Value.I")]
        colnames(result_mlm_spearman_A)[2]=paste0(vector[i],".mlm_spearman")
        colnames(result_mlm_spearman_I)[2]=paste0(vector[i],".mlm_spearman")
        rm(result_mlm_spearman)
      }
      if("mlm_lm" %in% enrichment_method){
        print("mlm_lm will start.")
        result_mlm_lm=lm_test(geneset=geneset,value=value)
        result_mlm_lm$two=TRUE
        result_mlm_lm$invert_A=ifelse(result_mlm_lm$r<=0,TRUE,FALSE)
        result_mlm_lm$invert_I=ifelse(result_mlm_lm$r>0,TRUE,FALSE)
        result_mlm_lm$P.Value.A=two2one(result_mlm_lm$P.Value, two = as.logical(result_mlm_lm$two), invert = as.logical(result_mlm_lm$invert_A))
        result_mlm_lm$P.Value.I=two2one(result_mlm_lm$P.Value, two = as.logical(result_mlm_lm$two), invert = as.logical(result_mlm_lm$invert_I))
        result_mlm_lm_A=result_mlm_lm[,c("geneset_id","P.Value.A")]
        result_mlm_lm_I=result_mlm_lm[,c("geneset_id","P.Value.I")]
        colnames(result_mlm_lm_A)[2]=paste0(vector[i],".mlm_lm")
        colnames(result_mlm_lm_I)[2]=paste0(vector[i],".mlm_lm")
        rm(result_mlm_lm)
      }
      if("mlm_biweight" %in% enrichment_method){
        print("mlm_biweight will start.")
        result_mlm_biweight=biweight_test(geneset=geneset,value=value)
        result_mlm_biweight$two=TRUE
        result_mlm_biweight$invert_A=ifelse(result_mlm_biweight$r<=0,TRUE,FALSE)
        result_mlm_biweight$invert_I=ifelse(result_mlm_biweight$r>0,TRUE,FALSE)
        result_mlm_biweight$P.Value.A=two2one(result_mlm_biweight$P.Value, two = as.logical(result_mlm_biweight$two), invert = as.logical(result_mlm_biweight$invert_A))
        result_mlm_biweight$P.Value.I=two2one(result_mlm_biweight$P.Value, two = as.logical(result_mlm_biweight$two), invert = as.logical(result_mlm_biweight$invert_I))
        result_mlm_biweight_A=result_mlm_biweight[,c("geneset_id","P.Value.A")]
        result_mlm_biweight_I=result_mlm_biweight[,c("geneset_id","P.Value.I")]
        colnames(result_mlm_biweight_A)[2]=paste0(vector[i],".mlm_biweight")
        colnames(result_mlm_biweight_I)[2]=paste0(vector[i],".mlm_biweight")
        rm(result_mlm_biweight)
      }
      if("mlm_distance" %in% enrichment_method){
        print("mlm_distance will start.")
        result_mlm_distance=distance_test(geneset=geneset,value=value)
        result_mlm_distance$two=TRUE
        result_mlm_distance$invert_A=ifelse(result_mlm_distance$r<=0,TRUE,FALSE)
        result_mlm_distance$invert_I=ifelse(result_mlm_distance$r>0,TRUE,FALSE)
        result_mlm_distance$P.Value.A=two2one(result_mlm_distance$P.Value, two = as.logical(result_mlm_distance$two), invert = as.logical(result_mlm_distance$invert_A))
        result_mlm_distance$P.Value.I=two2one(result_mlm_distance$P.Value, two = as.logical(result_mlm_distance$two), invert = as.logical(result_mlm_distance$invert_I))
        result_mlm_distance_A=result_mlm_distance[,c("geneset_id","P.Value.A")]
        result_mlm_distance_I=result_mlm_distance[,c("geneset_id","P.Value.I")]
        colnames(result_mlm_distance_A)[2]=paste0(vector[i],".mlm_distance")
        colnames(result_mlm_distance_I)[2]=paste0(vector[i],".mlm_distance")
        rm(result_mlm_distance)
      }
      if("mlm_percentage" %in% enrichment_method){
        print("mlm_percentage will start.")
        result_mlm_percentage=percentage_test(geneset=geneset,value=value)
        result_mlm_percentage$two=TRUE
        result_mlm_percentage$invert_A=ifelse(result_mlm_percentage$r<=0,TRUE,FALSE)
        result_mlm_percentage$invert_I=ifelse(result_mlm_percentage$r>0,TRUE,FALSE)
        result_mlm_percentage$P.Value.A=two2one(result_mlm_percentage$P.Value, two = as.logical(result_mlm_percentage$two), invert = as.logical(result_mlm_percentage$invert_A))
        result_mlm_percentage$P.Value.I=two2one(result_mlm_percentage$P.Value, two = as.logical(result_mlm_percentage$two), invert = as.logical(result_mlm_percentage$invert_I))
        result_mlm_percentage_A=result_mlm_percentage[,c("geneset_id","P.Value.A")]
        result_mlm_percentage_I=result_mlm_percentage[,c("geneset_id","P.Value.I")]
        colnames(result_mlm_percentage_A)[2]=paste0(vector[i],".mlm_percentage")
        colnames(result_mlm_percentage_I)[2]=paste0(vector[i],".mlm_percentage")
        rm(result_mlm_percentage)
      }
      if("mlm_blomqvist" %in% enrichment_method){
        print("mlm_blomqvist will start.")
        result_mlm_blomqvist=blomqvist_test(geneset=geneset,value=value)
        result_mlm_blomqvist$two=TRUE
        result_mlm_blomqvist$invert_A=ifelse(result_mlm_blomqvist$r<=0,TRUE,FALSE)
        result_mlm_blomqvist$invert_I=ifelse(result_mlm_blomqvist$r>0,TRUE,FALSE)
        result_mlm_blomqvist$P.Value.A=two2one(result_mlm_blomqvist$P.Value, two = as.logical(result_mlm_blomqvist$two), invert = as.logical(result_mlm_blomqvist$invert_A))
        result_mlm_blomqvist$P.Value.I=two2one(result_mlm_blomqvist$P.Value, two = as.logical(result_mlm_blomqvist$two), invert = as.logical(result_mlm_blomqvist$invert_I))
        result_mlm_blomqvist_A=result_mlm_blomqvist[,c("geneset_id","P.Value.A")]
        result_mlm_blomqvist_I=result_mlm_blomqvist[,c("geneset_id","P.Value.I")]
        colnames(result_mlm_blomqvist_A)[2]=paste0(vector[i],".mlm_blomqvist")
        colnames(result_mlm_blomqvist_I)[2]=paste0(vector[i],".mlm_blomqvist")
        rm(result_mlm_blomqvist)
      }
      if("mlm_hoeffding" %in% enrichment_method){
        print("mlm_hoeffding will start.")
        result_mlm_hoeffding=hoeffding_test(geneset=geneset,value=value)
        result_mlm_hoeffding$two=TRUE
        result_mlm_hoeffding$invert_A=ifelse(result_mlm_hoeffding$r<=0,TRUE,FALSE)
        result_mlm_hoeffding$invert_I=ifelse(result_mlm_hoeffding$r>0,TRUE,FALSE)
        result_mlm_hoeffding$P.Value.A=two2one(result_mlm_hoeffding$P.Value, two = as.logical(result_mlm_hoeffding$two), invert = as.logical(result_mlm_hoeffding$invert_A))
        result_mlm_hoeffding$P.Value.I=two2one(result_mlm_hoeffding$P.Value, two = as.logical(result_mlm_hoeffding$two), invert = as.logical(result_mlm_hoeffding$invert_I))
        result_mlm_hoeffding_A=result_mlm_hoeffding[,c("geneset_id","P.Value.A")]
        result_mlm_hoeffding_I=result_mlm_hoeffding[,c("geneset_id","P.Value.I")]
        colnames(result_mlm_hoeffding_A)[2]=paste0(vector[i],".mlm_hoeffding")
        colnames(result_mlm_hoeffding_I)[2]=paste0(vector[i],".mlm_hoeffding")
        rm(result_mlm_hoeffding)
      }
      if("mlm_gamma" %in% enrichment_method){
        print("mlm_gamma will start.")
        result_mlm_gamma=gamma_test(geneset=geneset,value=value)
        result_mlm_gamma$two=TRUE
        result_mlm_gamma$invert_A=ifelse(result_mlm_gamma$r<=0,TRUE,FALSE)
        result_mlm_gamma$invert_I=ifelse(result_mlm_gamma$r>0,TRUE,FALSE)
        result_mlm_gamma$P.Value.A=two2one(result_mlm_gamma$P.Value, two = as.logical(result_mlm_gamma$two), invert = as.logical(result_mlm_gamma$invert_A))
        result_mlm_gamma$P.Value.I=two2one(result_mlm_gamma$P.Value, two = as.logical(result_mlm_gamma$two), invert = as.logical(result_mlm_gamma$invert_I))
        result_mlm_gamma_A=result_mlm_gamma[,c("geneset_id","P.Value.A")]
        result_mlm_gamma_I=result_mlm_gamma[,c("geneset_id","P.Value.I")]
        colnames(result_mlm_gamma_A)[2]=paste0(vector[i],".mlm_gamma")
        colnames(result_mlm_gamma_I)[2]=paste0(vector[i],".mlm_gamma")
        rm(result_mlm_gamma)
      }
    }else {
      print("Not using mlm.")
    }
    ##########################################udt
    methods_to_check=c("udt_pearson","udt_kendall","udt_spearman","udt_lm","udt_biweight","udt_distance","udt_percentage","udt_blomqvist","udt_hoeffding","udt_gamma")
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
      if("udt_pearson" %in% enrichment_method){
        print("udt_pearson will start.")
        result_udt_pearson=pearson_test(geneset=geneset,value=value)
        result_udt_pearson$two=TRUE
        result_udt_pearson$invert_A=ifelse(result_udt_pearson$r<=0,TRUE,FALSE)
        result_udt_pearson$invert_I=ifelse(result_udt_pearson$r>0,TRUE,FALSE)
        result_udt_pearson$P.Value.A=two2one(result_udt_pearson$P.Value, two = as.logical(result_udt_pearson$two), invert = as.logical(result_udt_pearson$invert_A))
        result_udt_pearson$P.Value.I=two2one(result_udt_pearson$P.Value, two = as.logical(result_udt_pearson$two), invert = as.logical(result_udt_pearson$invert_I))
        result_udt_pearson_A=result_udt_pearson[,c("geneset_id","P.Value.A")]
        result_udt_pearson_I=result_udt_pearson[,c("geneset_id","P.Value.I")]
        colnames(result_udt_pearson_A)[2]=paste0(vector[i],".udt_pearson")
        colnames(result_udt_pearson_I)[2]=paste0(vector[i],".udt_pearson")
        rm(result_udt_pearson)
      }
      if("udt_kendall" %in% enrichment_method){
        print("udt_kendall will start.")
        result_udt_kendall=kendall_test(geneset=geneset,value=value)
        result_udt_kendall$two=TRUE
        result_udt_kendall$invert_A=ifelse(result_udt_kendall$r<=0,TRUE,FALSE)
        result_udt_kendall$invert_I=ifelse(result_udt_kendall$r>0,TRUE,FALSE)
        result_udt_kendall$P.Value.A=two2one(result_udt_kendall$P.Value, two = as.logical(result_udt_kendall$two), invert = as.logical(result_udt_kendall$invert_A))
        result_udt_kendall$P.Value.I=two2one(result_udt_kendall$P.Value, two = as.logical(result_udt_kendall$two), invert = as.logical(result_udt_kendall$invert_I))
        result_udt_kendall_A=result_udt_kendall[,c("geneset_id","P.Value.A")]
        result_udt_kendall_I=result_udt_kendall[,c("geneset_id","P.Value.I")]
        colnames(result_udt_kendall_A)[2]=paste0(vector[i],".udt_kendall")
        colnames(result_udt_kendall_I)[2]=paste0(vector[i],".udt_kendall")
        rm(result_udt_kendall)
      }
      if("udt_spearman" %in% enrichment_method){
        print("udt_spearman will start.")
        result_udt_spearman=spearman_test(geneset=geneset,value=value)
        result_udt_spearman$two=TRUE
        result_udt_spearman$invert_A=ifelse(result_udt_spearman$r<=0,TRUE,FALSE)
        result_udt_spearman$invert_I=ifelse(result_udt_spearman$r>0,TRUE,FALSE)
        result_udt_spearman$P.Value.A=two2one(result_udt_spearman$P.Value, two = as.logical(result_udt_spearman$two), invert = as.logical(result_udt_spearman$invert_A))
        result_udt_spearman$P.Value.I=two2one(result_udt_spearman$P.Value, two = as.logical(result_udt_spearman$two), invert = as.logical(result_udt_spearman$invert_I))
        result_udt_spearman_A=result_udt_spearman[,c("geneset_id","P.Value.A")]
        result_udt_spearman_I=result_udt_spearman[,c("geneset_id","P.Value.I")]
        colnames(result_udt_spearman_A)[2]=paste0(vector[i],".udt_spearman")
        colnames(result_udt_spearman_I)[2]=paste0(vector[i],".udt_spearman")
        rm(result_udt_spearman)
      }
      if("udt_lm" %in% enrichment_method){
        print("udt_lm will start.")
        result_udt_lm=lm_test(geneset=geneset,value=value)
        result_udt_lm$two=TRUE
        result_udt_lm$invert_A=ifelse(result_udt_lm$r<=0,TRUE,FALSE)
        result_udt_lm$invert_I=ifelse(result_udt_lm$r>0,TRUE,FALSE)
        result_udt_lm$P.Value.A=two2one(result_udt_lm$P.Value, two = as.logical(result_udt_lm$two), invert = as.logical(result_udt_lm$invert_A))
        result_udt_lm$P.Value.I=two2one(result_udt_lm$P.Value, two = as.logical(result_udt_lm$two), invert = as.logical(result_udt_lm$invert_I))
        result_udt_lm_A=result_udt_lm[,c("geneset_id","P.Value.A")]
        result_udt_lm_I=result_udt_lm[,c("geneset_id","P.Value.I")]
        colnames(result_udt_lm_A)[2]=paste0(vector[i],".udt_lm")
        colnames(result_udt_lm_I)[2]=paste0(vector[i],".udt_lm")
        rm(result_udt_lm)
      }
      if("udt_biweight" %in% enrichment_method){
        print("udt_biweight will start.")
        result_udt_biweight=biweight_test(geneset=geneset,value=value)
        result_udt_biweight$two=TRUE
        result_udt_biweight$invert_A=ifelse(result_udt_biweight$r<=0,TRUE,FALSE)
        result_udt_biweight$invert_I=ifelse(result_udt_biweight$r>0,TRUE,FALSE)
        result_udt_biweight$P.Value.A=two2one(result_udt_biweight$P.Value, two = as.logical(result_udt_biweight$two), invert = as.logical(result_udt_biweight$invert_A))
        result_udt_biweight$P.Value.I=two2one(result_udt_biweight$P.Value, two = as.logical(result_udt_biweight$two), invert = as.logical(result_udt_biweight$invert_I))
        result_udt_biweight_A=result_udt_biweight[,c("geneset_id","P.Value.A")]
        result_udt_biweight_I=result_udt_biweight[,c("geneset_id","P.Value.I")]
        colnames(result_udt_biweight_A)[2]=paste0(vector[i],".udt_biweight")
        colnames(result_udt_biweight_I)[2]=paste0(vector[i],".udt_biweight")
        rm(result_udt_biweight)
      }
      if("udt_distance" %in% enrichment_method){
        print("udt_distance will start.")
        result_udt_distance=distance_test(geneset=geneset,value=value)
        result_udt_distance$two=TRUE
        result_udt_distance$invert_A=ifelse(result_udt_distance$r<=0,TRUE,FALSE)
        result_udt_distance$invert_I=ifelse(result_udt_distance$r>0,TRUE,FALSE)
        result_udt_distance$P.Value.A=two2one(result_udt_distance$P.Value, two = as.logical(result_udt_distance$two), invert = as.logical(result_udt_distance$invert_A))
        result_udt_distance$P.Value.I=two2one(result_udt_distance$P.Value, two = as.logical(result_udt_distance$two), invert = as.logical(result_udt_distance$invert_I))
        result_udt_distance_A=result_udt_distance[,c("geneset_id","P.Value.A")]
        result_udt_distance_I=result_udt_distance[,c("geneset_id","P.Value.I")]
        colnames(result_udt_distance_A)[2]=paste0(vector[i],".udt_distance")
        colnames(result_udt_distance_I)[2]=paste0(vector[i],".udt_distance")
        rm(result_udt_distance)
      }
      if("udt_percentage" %in% enrichment_method){
        print("udt_percentage will start.")
        result_udt_percentage=percentage_test(geneset=geneset,value=value)
        result_udt_percentage$two=TRUE
        result_udt_percentage$invert_A=ifelse(result_udt_percentage$r<=0,TRUE,FALSE)
        result_udt_percentage$invert_I=ifelse(result_udt_percentage$r>0,TRUE,FALSE)
        result_udt_percentage$P.Value.A=two2one(result_udt_percentage$P.Value, two = as.logical(result_udt_percentage$two), invert = as.logical(result_udt_percentage$invert_A))
        result_udt_percentage$P.Value.I=two2one(result_udt_percentage$P.Value, two = as.logical(result_udt_percentage$two), invert = as.logical(result_udt_percentage$invert_I))
        result_udt_percentage_A=result_udt_percentage[,c("geneset_id","P.Value.A")]
        result_udt_percentage_I=result_udt_percentage[,c("geneset_id","P.Value.I")]
        colnames(result_udt_percentage_A)[2]=paste0(vector[i],".udt_percentage")
        colnames(result_udt_percentage_I)[2]=paste0(vector[i],".udt_percentage")
        rm(result_udt_percentage)
      }
      if("udt_blomqvist" %in% enrichment_method){
        print("udt_blomqvist will start.")
        result_udt_blomqvist=blomqvist_test(geneset=geneset,value=value)
        result_udt_blomqvist$two=TRUE
        result_udt_blomqvist$invert_A=ifelse(result_udt_blomqvist$r<=0,TRUE,FALSE)
        result_udt_blomqvist$invert_I=ifelse(result_udt_blomqvist$r>0,TRUE,FALSE)
        result_udt_blomqvist$P.Value.A=two2one(result_udt_blomqvist$P.Value, two = as.logical(result_udt_blomqvist$two), invert = as.logical(result_udt_blomqvist$invert_A))
        result_udt_blomqvist$P.Value.I=two2one(result_udt_blomqvist$P.Value, two = as.logical(result_udt_blomqvist$two), invert = as.logical(result_udt_blomqvist$invert_I))
        result_udt_blomqvist_A=result_udt_blomqvist[,c("geneset_id","P.Value.A")]
        result_udt_blomqvist_I=result_udt_blomqvist[,c("geneset_id","P.Value.I")]
        colnames(result_udt_blomqvist_A)[2]=paste0(vector[i],".udt_blomqvist")
        colnames(result_udt_blomqvist_I)[2]=paste0(vector[i],".udt_blomqvist")
        rm(result_udt_blomqvist)
      }
      if("udt_hoeffding" %in% enrichment_method){
        print("udt_hoeffding will start.")
        result_udt_hoeffding=hoeffding_test(geneset=geneset,value=value)
        result_udt_hoeffding$two=TRUE
        result_udt_hoeffding$invert_A=ifelse(result_udt_hoeffding$r<=0,TRUE,FALSE)
        result_udt_hoeffding$invert_I=ifelse(result_udt_hoeffding$r>0,TRUE,FALSE)
        result_udt_hoeffding$P.Value.A=two2one(result_udt_hoeffding$P.Value, two = as.logical(result_udt_hoeffding$two), invert = as.logical(result_udt_hoeffding$invert_A))
        result_udt_hoeffding$P.Value.I=two2one(result_udt_hoeffding$P.Value, two = as.logical(result_udt_hoeffding$two), invert = as.logical(result_udt_hoeffding$invert_I))
        result_udt_hoeffding_A=result_udt_hoeffding[,c("geneset_id","P.Value.A")]
        result_udt_hoeffding_I=result_udt_hoeffding[,c("geneset_id","P.Value.I")]
        colnames(result_udt_hoeffding_A)[2]=paste0(vector[i],".udt_hoeffding")
        colnames(result_udt_hoeffding_I)[2]=paste0(vector[i],".udt_hoeffding")
        rm(result_udt_hoeffding)
      }
      if("udt_gamma" %in% enrichment_method){
        print("udt_gamma will start.")
        result_udt_gamma=gamma_test(geneset=geneset,value=value)
        result_udt_gamma$two=TRUE
        result_udt_gamma$invert_A=ifelse(result_udt_gamma$r<=0,TRUE,FALSE)
        result_udt_gamma$invert_I=ifelse(result_udt_gamma$r>0,TRUE,FALSE)
        result_udt_gamma$P.Value.A=two2one(result_udt_gamma$P.Value, two = as.logical(result_udt_gamma$two), invert = as.logical(result_udt_gamma$invert_A))
        result_udt_gamma$P.Value.I=two2one(result_udt_gamma$P.Value, two = as.logical(result_udt_gamma$two), invert = as.logical(result_udt_gamma$invert_I))
        result_udt_gamma_A=result_udt_gamma[,c("geneset_id","P.Value.A")]
        result_udt_gamma_I=result_udt_gamma[,c("geneset_id","P.Value.I")]
        colnames(result_udt_gamma_A)[2]=paste0(vector[i],".udt_gamma")
        colnames(result_udt_gamma_I)[2]=paste0(vector[i],".udt_gamma")
        rm(result_udt_gamma)
      }
    }else {
      print("Not using udt.")
    }
    ##########################################ulm
    methods_to_check=c("ulm_pearson","ulm_kendall","ulm_spearman","ulm_lm","ulm_biweight","ulm_distance","ulm_percentage","ulm_blomqvist","ulm_hoeffding","ulm_gamma")
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
      if("ulm_pearson" %in% enrichment_method){
        print("ulm_pearson will start.")
        result_ulm_pearson=pearson_test(geneset=geneset,value=value)
        result_ulm_pearson$two=TRUE
        result_ulm_pearson$invert_A=ifelse(result_ulm_pearson$r<=0,TRUE,FALSE)
        result_ulm_pearson$invert_I=ifelse(result_ulm_pearson$r>0,TRUE,FALSE)
        result_ulm_pearson$P.Value.A=two2one(result_ulm_pearson$P.Value, two = as.logical(result_ulm_pearson$two), invert = as.logical(result_ulm_pearson$invert_A))
        result_ulm_pearson$P.Value.I=two2one(result_ulm_pearson$P.Value, two = as.logical(result_ulm_pearson$two), invert = as.logical(result_ulm_pearson$invert_I))
        result_ulm_pearson_A=result_ulm_pearson[,c("geneset_id","P.Value.A")]
        result_ulm_pearson_I=result_ulm_pearson[,c("geneset_id","P.Value.I")]
        colnames(result_ulm_pearson_A)[2]=paste0(vector[i],".ulm_pearson")
        colnames(result_ulm_pearson_I)[2]=paste0(vector[i],".ulm_pearson")
        rm(result_ulm_pearson)
      }
      if("ulm_kendall" %in% enrichment_method){
        print("ulm_kendall will start.")
        result_ulm_kendall=kendall_test(geneset=geneset,value=value)
        result_ulm_kendall$two=TRUE
        result_ulm_kendall$invert_A=ifelse(result_ulm_kendall$r<=0,TRUE,FALSE)
        result_ulm_kendall$invert_I=ifelse(result_ulm_kendall$r>0,TRUE,FALSE)
        result_ulm_kendall$P.Value.A=two2one(result_ulm_kendall$P.Value, two = as.logical(result_ulm_kendall$two), invert = as.logical(result_ulm_kendall$invert_A))
        result_ulm_kendall$P.Value.I=two2one(result_ulm_kendall$P.Value, two = as.logical(result_ulm_kendall$two), invert = as.logical(result_ulm_kendall$invert_I))
        result_ulm_kendall_A=result_ulm_kendall[,c("geneset_id","P.Value.A")]
        result_ulm_kendall_I=result_ulm_kendall[,c("geneset_id","P.Value.I")]
        colnames(result_ulm_kendall_A)[2]=paste0(vector[i],".ulm_kendall")
        colnames(result_ulm_kendall_I)[2]=paste0(vector[i],".ulm_kendall")
        rm(result_ulm_kendall)
      }
      if("ulm_spearman" %in% enrichment_method){
        print("ulm_spearman will start.")
        result_ulm_spearman=spearman_test(geneset=geneset,value=value)
        result_ulm_spearman$two=TRUE
        result_ulm_spearman$invert_A=ifelse(result_ulm_spearman$r<=0,TRUE,FALSE)
        result_ulm_spearman$invert_I=ifelse(result_ulm_spearman$r>0,TRUE,FALSE)
        result_ulm_spearman$P.Value.A=two2one(result_ulm_spearman$P.Value, two = as.logical(result_ulm_spearman$two), invert = as.logical(result_ulm_spearman$invert_A))
        result_ulm_spearman$P.Value.I=two2one(result_ulm_spearman$P.Value, two = as.logical(result_ulm_spearman$two), invert = as.logical(result_ulm_spearman$invert_I))
        result_ulm_spearman_A=result_ulm_spearman[,c("geneset_id","P.Value.A")]
        result_ulm_spearman_I=result_ulm_spearman[,c("geneset_id","P.Value.I")]
        colnames(result_ulm_spearman_A)[2]=paste0(vector[i],".ulm_spearman")
        colnames(result_ulm_spearman_I)[2]=paste0(vector[i],".ulm_spearman")
        rm(result_ulm_spearman)
      }
      if("ulm_lm" %in% enrichment_method){
        print("ulm_lm will start.")
        result_ulm_lm=lm_test(geneset=geneset,value=value)
        result_ulm_lm$two=TRUE
        result_ulm_lm$invert_A=ifelse(result_ulm_lm$r<=0,TRUE,FALSE)
        result_ulm_lm$invert_I=ifelse(result_ulm_lm$r>0,TRUE,FALSE)
        result_ulm_lm$P.Value.A=two2one(result_ulm_lm$P.Value, two = as.logical(result_ulm_lm$two), invert = as.logical(result_ulm_lm$invert_A))
        result_ulm_lm$P.Value.I=two2one(result_ulm_lm$P.Value, two = as.logical(result_ulm_lm$two), invert = as.logical(result_ulm_lm$invert_I))
        result_ulm_lm_A=result_ulm_lm[,c("geneset_id","P.Value.A")]
        result_ulm_lm_I=result_ulm_lm[,c("geneset_id","P.Value.I")]
        colnames(result_ulm_lm_A)[2]=paste0(vector[i],".ulm_lm")
        colnames(result_ulm_lm_I)[2]=paste0(vector[i],".ulm_lm")
        rm(result_ulm_lm)
      }
      if("ulm_biweight" %in% enrichment_method){
        print("ulm_biweight will start.")
        result_ulm_biweight=biweight_test(geneset=geneset,value=value)
        result_ulm_biweight$two=TRUE
        result_ulm_biweight$invert_A=ifelse(result_ulm_biweight$r<=0,TRUE,FALSE)
        result_ulm_biweight$invert_I=ifelse(result_ulm_biweight$r>0,TRUE,FALSE)
        result_ulm_biweight$P.Value.A=two2one(result_ulm_biweight$P.Value, two = as.logical(result_ulm_biweight$two), invert = as.logical(result_ulm_biweight$invert_A))
        result_ulm_biweight$P.Value.I=two2one(result_ulm_biweight$P.Value, two = as.logical(result_ulm_biweight$two), invert = as.logical(result_ulm_biweight$invert_I))
        result_ulm_biweight_A=result_ulm_biweight[,c("geneset_id","P.Value.A")]
        result_ulm_biweight_I=result_ulm_biweight[,c("geneset_id","P.Value.I")]
        colnames(result_ulm_biweight_A)[2]=paste0(vector[i],".ulm_biweight")
        colnames(result_ulm_biweight_I)[2]=paste0(vector[i],".ulm_biweight")
        rm(result_ulm_biweight)
      }
      if("ulm_distance" %in% enrichment_method){
        print("ulm_distance will start.")
        result_ulm_distance=distance_test(geneset=geneset,value=value)
        result_ulm_distance$two=TRUE
        result_ulm_distance$invert_A=ifelse(result_ulm_distance$r<=0,TRUE,FALSE)
        result_ulm_distance$invert_I=ifelse(result_ulm_distance$r>0,TRUE,FALSE)
        result_ulm_distance$P.Value.A=two2one(result_ulm_distance$P.Value, two = as.logical(result_ulm_distance$two), invert = as.logical(result_ulm_distance$invert_A))
        result_ulm_distance$P.Value.I=two2one(result_ulm_distance$P.Value, two = as.logical(result_ulm_distance$two), invert = as.logical(result_ulm_distance$invert_I))
        result_ulm_distance_A=result_ulm_distance[,c("geneset_id","P.Value.A")]
        result_ulm_distance_I=result_ulm_distance[,c("geneset_id","P.Value.I")]
        colnames(result_ulm_distance_A)[2]=paste0(vector[i],".ulm_distance")
        colnames(result_ulm_distance_I)[2]=paste0(vector[i],".ulm_distance")
        rm(result_ulm_distance)
      }
      if("ulm_percentage" %in% enrichment_method){
        print("ulm_percentage will start.")
        result_ulm_percentage=percentage_test(geneset=geneset,value=value)
        result_ulm_percentage$two=TRUE
        result_ulm_percentage$invert_A=ifelse(result_ulm_percentage$r<=0,TRUE,FALSE)
        result_ulm_percentage$invert_I=ifelse(result_ulm_percentage$r>0,TRUE,FALSE)
        result_ulm_percentage$P.Value.A=two2one(result_ulm_percentage$P.Value, two = as.logical(result_ulm_percentage$two), invert = as.logical(result_ulm_percentage$invert_A))
        result_ulm_percentage$P.Value.I=two2one(result_ulm_percentage$P.Value, two = as.logical(result_ulm_percentage$two), invert = as.logical(result_ulm_percentage$invert_I))
        result_ulm_percentage_A=result_ulm_percentage[,c("geneset_id","P.Value.A")]
        result_ulm_percentage_I=result_ulm_percentage[,c("geneset_id","P.Value.I")]
        colnames(result_ulm_percentage_A)[2]=paste0(vector[i],".ulm_percentage")
        colnames(result_ulm_percentage_I)[2]=paste0(vector[i],".ulm_percentage")
        rm(result_ulm_percentage)
      }
      if("ulm_blomqvist" %in% enrichment_method){
        print("ulm_blomqvist will start.")
        result_ulm_blomqvist=blomqvist_test(geneset=geneset,value=value)
        result_ulm_blomqvist$two=TRUE
        result_ulm_blomqvist$invert_A=ifelse(result_ulm_blomqvist$r<=0,TRUE,FALSE)
        result_ulm_blomqvist$invert_I=ifelse(result_ulm_blomqvist$r>0,TRUE,FALSE)
        result_ulm_blomqvist$P.Value.A=two2one(result_ulm_blomqvist$P.Value, two = as.logical(result_ulm_blomqvist$two), invert = as.logical(result_ulm_blomqvist$invert_A))
        result_ulm_blomqvist$P.Value.I=two2one(result_ulm_blomqvist$P.Value, two = as.logical(result_ulm_blomqvist$two), invert = as.logical(result_ulm_blomqvist$invert_I))
        result_ulm_blomqvist_A=result_ulm_blomqvist[,c("geneset_id","P.Value.A")]
        result_ulm_blomqvist_I=result_ulm_blomqvist[,c("geneset_id","P.Value.I")]
        colnames(result_ulm_blomqvist_A)[2]=paste0(vector[i],".ulm_blomqvist")
        colnames(result_ulm_blomqvist_I)[2]=paste0(vector[i],".ulm_blomqvist")
        rm(result_ulm_blomqvist)
      }
      if("ulm_hoeffding" %in% enrichment_method){
        print("ulm_hoeffding will start.")
        result_ulm_hoeffding=hoeffding_test(geneset=geneset,value=value)
        result_ulm_hoeffding$two=TRUE
        result_ulm_hoeffding$invert_A=ifelse(result_ulm_hoeffding$r<=0,TRUE,FALSE)
        result_ulm_hoeffding$invert_I=ifelse(result_ulm_hoeffding$r>0,TRUE,FALSE)
        result_ulm_hoeffding$P.Value.A=two2one(result_ulm_hoeffding$P.Value, two = as.logical(result_ulm_hoeffding$two), invert = as.logical(result_ulm_hoeffding$invert_A))
        result_ulm_hoeffding$P.Value.I=two2one(result_ulm_hoeffding$P.Value, two = as.logical(result_ulm_hoeffding$two), invert = as.logical(result_ulm_hoeffding$invert_I))
        result_ulm_hoeffding_A=result_ulm_hoeffding[,c("geneset_id","P.Value.A")]
        result_ulm_hoeffding_I=result_ulm_hoeffding[,c("geneset_id","P.Value.I")]
        colnames(result_ulm_hoeffding_A)[2]=paste0(vector[i],".ulm_hoeffding")
        colnames(result_ulm_hoeffding_I)[2]=paste0(vector[i],".ulm_hoeffding")
        rm(result_ulm_hoeffding)
      }
      if("ulm_gamma" %in% enrichment_method){
        print("ulm_gamma will start.")
        result_ulm_gamma=gamma_test(geneset=geneset,value=value)
        result_ulm_gamma$two=TRUE
        result_ulm_gamma$invert_A=ifelse(result_ulm_gamma$r<=0,TRUE,FALSE)
        result_ulm_gamma$invert_I=ifelse(result_ulm_gamma$r>0,TRUE,FALSE)
        result_ulm_gamma$P.Value.A=two2one(result_ulm_gamma$P.Value, two = as.logical(result_ulm_gamma$two), invert = as.logical(result_ulm_gamma$invert_A))
        result_ulm_gamma$P.Value.I=two2one(result_ulm_gamma$P.Value, two = as.logical(result_ulm_gamma$two), invert = as.logical(result_ulm_gamma$invert_I))
        result_ulm_gamma_A=result_ulm_gamma[,c("geneset_id","P.Value.A")]
        result_ulm_gamma_I=result_ulm_gamma[,c("geneset_id","P.Value.I")]
        colnames(result_ulm_gamma_A)[2]=paste0(vector[i],".ulm_gamma")
        colnames(result_ulm_gamma_I)[2]=paste0(vector[i],".ulm_gamma")
        rm(result_ulm_gamma)
      }
    }else {
      print("Not using ulm.")
    }
    ##########################################viper
    methods_to_check=c("viper_pearson","viper_kendall","viper_spearman","viper_lm","viper_biweight","viper_distance","viper_percentage","viper_blomqvist","viper_hoeffding","viper_gamma")
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
      if("viper_pearson" %in% enrichment_method){
        print("viper_pearson will start.")
        result_viper_pearson=pearson_test(geneset=geneset,value=value)
        result_viper_pearson$two=TRUE
        result_viper_pearson$invert_A=ifelse(result_viper_pearson$r<=0,TRUE,FALSE)
        result_viper_pearson$invert_I=ifelse(result_viper_pearson$r>0,TRUE,FALSE)
        result_viper_pearson$P.Value.A=two2one(result_viper_pearson$P.Value, two = as.logical(result_viper_pearson$two), invert = as.logical(result_viper_pearson$invert_A))
        result_viper_pearson$P.Value.I=two2one(result_viper_pearson$P.Value, two = as.logical(result_viper_pearson$two), invert = as.logical(result_viper_pearson$invert_I))
        result_viper_pearson_A=result_viper_pearson[,c("geneset_id","P.Value.A")]
        result_viper_pearson_I=result_viper_pearson[,c("geneset_id","P.Value.I")]
        colnames(result_viper_pearson_A)[2]=paste0(vector[i],".viper_pearson")
        colnames(result_viper_pearson_I)[2]=paste0(vector[i],".viper_pearson")
        rm(result_viper_pearson)
      }
      if("viper_kendall" %in% enrichment_method){
        print("viper_kendall will start.")
        result_viper_kendall=kendall_test(geneset=geneset,value=value)
        result_viper_kendall$two=TRUE
        result_viper_kendall$invert_A=ifelse(result_viper_kendall$r<=0,TRUE,FALSE)
        result_viper_kendall$invert_I=ifelse(result_viper_kendall$r>0,TRUE,FALSE)
        result_viper_kendall$P.Value.A=two2one(result_viper_kendall$P.Value, two = as.logical(result_viper_kendall$two), invert = as.logical(result_viper_kendall$invert_A))
        result_viper_kendall$P.Value.I=two2one(result_viper_kendall$P.Value, two = as.logical(result_viper_kendall$two), invert = as.logical(result_viper_kendall$invert_I))
        result_viper_kendall_A=result_viper_kendall[,c("geneset_id","P.Value.A")]
        result_viper_kendall_I=result_viper_kendall[,c("geneset_id","P.Value.I")]
        colnames(result_viper_kendall_A)[2]=paste0(vector[i],".viper_kendall")
        colnames(result_viper_kendall_I)[2]=paste0(vector[i],".viper_kendall")
        rm(result_viper_kendall)
      }
      if("viper_spearman" %in% enrichment_method){
        print("viper_spearman will start.")
        result_viper_spearman=spearman_test(geneset=geneset,value=value)
        result_viper_spearman$two=TRUE
        result_viper_spearman$invert_A=ifelse(result_viper_spearman$r<=0,TRUE,FALSE)
        result_viper_spearman$invert_I=ifelse(result_viper_spearman$r>0,TRUE,FALSE)
        result_viper_spearman$P.Value.A=two2one(result_viper_spearman$P.Value, two = as.logical(result_viper_spearman$two), invert = as.logical(result_viper_spearman$invert_A))
        result_viper_spearman$P.Value.I=two2one(result_viper_spearman$P.Value, two = as.logical(result_viper_spearman$two), invert = as.logical(result_viper_spearman$invert_I))
        result_viper_spearman_A=result_viper_spearman[,c("geneset_id","P.Value.A")]
        result_viper_spearman_I=result_viper_spearman[,c("geneset_id","P.Value.I")]
        colnames(result_viper_spearman_A)[2]=paste0(vector[i],".viper_spearman")
        colnames(result_viper_spearman_I)[2]=paste0(vector[i],".viper_spearman")
        rm(result_viper_spearman)
      }
      if("viper_lm" %in% enrichment_method){
        print("viper_lm will start.")
        result_viper_lm=lm_test(geneset=geneset,value=value)
        result_viper_lm$two=TRUE
        result_viper_lm$invert_A=ifelse(result_viper_lm$r<=0,TRUE,FALSE)
        result_viper_lm$invert_I=ifelse(result_viper_lm$r>0,TRUE,FALSE)
        result_viper_lm$P.Value.A=two2one(result_viper_lm$P.Value, two = as.logical(result_viper_lm$two), invert = as.logical(result_viper_lm$invert_A))
        result_viper_lm$P.Value.I=two2one(result_viper_lm$P.Value, two = as.logical(result_viper_lm$two), invert = as.logical(result_viper_lm$invert_I))
        result_viper_lm_A=result_viper_lm[,c("geneset_id","P.Value.A")]
        result_viper_lm_I=result_viper_lm[,c("geneset_id","P.Value.I")]
        colnames(result_viper_lm_A)[2]=paste0(vector[i],".viper_lm")
        colnames(result_viper_lm_I)[2]=paste0(vector[i],".viper_lm")
        rm(result_viper_lm)
      }
      if("viper_biweight" %in% enrichment_method){
        print("viper_biweight will start.")
        result_viper_biweight=biweight_test(geneset=geneset,value=value)
        result_viper_biweight$two=TRUE
        result_viper_biweight$invert_A=ifelse(result_viper_biweight$r<=0,TRUE,FALSE)
        result_viper_biweight$invert_I=ifelse(result_viper_biweight$r>0,TRUE,FALSE)
        result_viper_biweight$P.Value.A=two2one(result_viper_biweight$P.Value, two = as.logical(result_viper_biweight$two), invert = as.logical(result_viper_biweight$invert_A))
        result_viper_biweight$P.Value.I=two2one(result_viper_biweight$P.Value, two = as.logical(result_viper_biweight$two), invert = as.logical(result_viper_biweight$invert_I))
        result_viper_biweight_A=result_viper_biweight[,c("geneset_id","P.Value.A")]
        result_viper_biweight_I=result_viper_biweight[,c("geneset_id","P.Value.I")]
        colnames(result_viper_biweight_A)[2]=paste0(vector[i],".viper_biweight")
        colnames(result_viper_biweight_I)[2]=paste0(vector[i],".viper_biweight")
        rm(result_viper_biweight)
      }
      if("viper_distance" %in% enrichment_method){
        print("viper_distance will start.")
        result_viper_distance=distance_test(geneset=geneset,value=value)
        result_viper_distance$two=TRUE
        result_viper_distance$invert_A=ifelse(result_viper_distance$r<=0,TRUE,FALSE)
        result_viper_distance$invert_I=ifelse(result_viper_distance$r>0,TRUE,FALSE)
        result_viper_distance$P.Value.A=two2one(result_viper_distance$P.Value, two = as.logical(result_viper_distance$two), invert = as.logical(result_viper_distance$invert_A))
        result_viper_distance$P.Value.I=two2one(result_viper_distance$P.Value, two = as.logical(result_viper_distance$two), invert = as.logical(result_viper_distance$invert_I))
        result_viper_distance_A=result_viper_distance[,c("geneset_id","P.Value.A")]
        result_viper_distance_I=result_viper_distance[,c("geneset_id","P.Value.I")]
        colnames(result_viper_distance_A)[2]=paste0(vector[i],".viper_distance")
        colnames(result_viper_distance_I)[2]=paste0(vector[i],".viper_distance")
        rm(result_viper_distance)
      }
      if("viper_percentage" %in% enrichment_method){
        print("viper_percentage will start.")
        result_viper_percentage=percentage_test(geneset=geneset,value=value)
        result_viper_percentage$two=TRUE
        result_viper_percentage$invert_A=ifelse(result_viper_percentage$r<=0,TRUE,FALSE)
        result_viper_percentage$invert_I=ifelse(result_viper_percentage$r>0,TRUE,FALSE)
        result_viper_percentage$P.Value.A=two2one(result_viper_percentage$P.Value, two = as.logical(result_viper_percentage$two), invert = as.logical(result_viper_percentage$invert_A))
        result_viper_percentage$P.Value.I=two2one(result_viper_percentage$P.Value, two = as.logical(result_viper_percentage$two), invert = as.logical(result_viper_percentage$invert_I))
        result_viper_percentage_A=result_viper_percentage[,c("geneset_id","P.Value.A")]
        result_viper_percentage_I=result_viper_percentage[,c("geneset_id","P.Value.I")]
        colnames(result_viper_percentage_A)[2]=paste0(vector[i],".viper_percentage")
        colnames(result_viper_percentage_I)[2]=paste0(vector[i],".viper_percentage")
        rm(result_viper_percentage)
      }
      if("viper_blomqvist" %in% enrichment_method){
        print("viper_blomqvist will start.")
        result_viper_blomqvist=blomqvist_test(geneset=geneset,value=value)
        result_viper_blomqvist$two=TRUE
        result_viper_blomqvist$invert_A=ifelse(result_viper_blomqvist$r<=0,TRUE,FALSE)
        result_viper_blomqvist$invert_I=ifelse(result_viper_blomqvist$r>0,TRUE,FALSE)
        result_viper_blomqvist$P.Value.A=two2one(result_viper_blomqvist$P.Value, two = as.logical(result_viper_blomqvist$two), invert = as.logical(result_viper_blomqvist$invert_A))
        result_viper_blomqvist$P.Value.I=two2one(result_viper_blomqvist$P.Value, two = as.logical(result_viper_blomqvist$two), invert = as.logical(result_viper_blomqvist$invert_I))
        result_viper_blomqvist_A=result_viper_blomqvist[,c("geneset_id","P.Value.A")]
        result_viper_blomqvist_I=result_viper_blomqvist[,c("geneset_id","P.Value.I")]
        colnames(result_viper_blomqvist_A)[2]=paste0(vector[i],".viper_blomqvist")
        colnames(result_viper_blomqvist_I)[2]=paste0(vector[i],".viper_blomqvist")
        rm(result_viper_blomqvist)
      }
      if("viper_hoeffding" %in% enrichment_method){
        print("viper_hoeffding will start.")
        result_viper_hoeffding=hoeffding_test(geneset=geneset,value=value)
        result_viper_hoeffding$two=TRUE
        result_viper_hoeffding$invert_A=ifelse(result_viper_hoeffding$r<=0,TRUE,FALSE)
        result_viper_hoeffding$invert_I=ifelse(result_viper_hoeffding$r>0,TRUE,FALSE)
        result_viper_hoeffding$P.Value.A=two2one(result_viper_hoeffding$P.Value, two = as.logical(result_viper_hoeffding$two), invert = as.logical(result_viper_hoeffding$invert_A))
        result_viper_hoeffding$P.Value.I=two2one(result_viper_hoeffding$P.Value, two = as.logical(result_viper_hoeffding$two), invert = as.logical(result_viper_hoeffding$invert_I))
        result_viper_hoeffding_A=result_viper_hoeffding[,c("geneset_id","P.Value.A")]
        result_viper_hoeffding_I=result_viper_hoeffding[,c("geneset_id","P.Value.I")]
        colnames(result_viper_hoeffding_A)[2]=paste0(vector[i],".viper_hoeffding")
        colnames(result_viper_hoeffding_I)[2]=paste0(vector[i],".viper_hoeffding")
        rm(result_viper_hoeffding)
      }
      if("viper_gamma" %in% enrichment_method){
        print("viper_gamma will start.")
        result_viper_gamma=gamma_test(geneset=geneset,value=value)
        result_viper_gamma$two=TRUE
        result_viper_gamma$invert_A=ifelse(result_viper_gamma$r<=0,TRUE,FALSE)
        result_viper_gamma$invert_I=ifelse(result_viper_gamma$r>0,TRUE,FALSE)
        result_viper_gamma$P.Value.A=two2one(result_viper_gamma$P.Value, two = as.logical(result_viper_gamma$two), invert = as.logical(result_viper_gamma$invert_A))
        result_viper_gamma$P.Value.I=two2one(result_viper_gamma$P.Value, two = as.logical(result_viper_gamma$two), invert = as.logical(result_viper_gamma$invert_I))
        result_viper_gamma_A=result_viper_gamma[,c("geneset_id","P.Value.A")]
        result_viper_gamma_I=result_viper_gamma[,c("geneset_id","P.Value.I")]
        colnames(result_viper_gamma_A)[2]=paste0(vector[i],".viper_gamma")
        colnames(result_viper_gamma_I)[2]=paste0(vector[i],".viper_gamma")
        rm(result_viper_gamma)
      }
    }else {
      print("Not using viper.")
    }
    ##########################################wmean
    methods_to_check=c("wmean_pearson","wmean_kendall","wmean_spearman","wmean_lm","wmean_biweight","wmean_distance","wmean_percentage","wmean_blomqvist","wmean_hoeffding","wmean_gamma")
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
      if("wmean_pearson" %in% enrichment_method){
        print("wmean_pearson will start.")
        result_wmean_pearson=pearson_test(geneset=geneset,value=value)
        result_wmean_pearson$two=TRUE
        result_wmean_pearson$invert_A=ifelse(result_wmean_pearson$r<=0,TRUE,FALSE)
        result_wmean_pearson$invert_I=ifelse(result_wmean_pearson$r>0,TRUE,FALSE)
        result_wmean_pearson$P.Value.A=two2one(result_wmean_pearson$P.Value, two = as.logical(result_wmean_pearson$two), invert = as.logical(result_wmean_pearson$invert_A))
        result_wmean_pearson$P.Value.I=two2one(result_wmean_pearson$P.Value, two = as.logical(result_wmean_pearson$two), invert = as.logical(result_wmean_pearson$invert_I))
        result_wmean_pearson_A=result_wmean_pearson[,c("geneset_id","P.Value.A")]
        result_wmean_pearson_I=result_wmean_pearson[,c("geneset_id","P.Value.I")]
        colnames(result_wmean_pearson_A)[2]=paste0(vector[i],".wmean_pearson")
        colnames(result_wmean_pearson_I)[2]=paste0(vector[i],".wmean_pearson")
        rm(result_wmean_pearson)
      }
      if("wmean_kendall" %in% enrichment_method){
        print("wmean_kendall will start.")
        result_wmean_kendall=kendall_test(geneset=geneset,value=value)
        result_wmean_kendall$two=TRUE
        result_wmean_kendall$invert_A=ifelse(result_wmean_kendall$r<=0,TRUE,FALSE)
        result_wmean_kendall$invert_I=ifelse(result_wmean_kendall$r>0,TRUE,FALSE)
        result_wmean_kendall$P.Value.A=two2one(result_wmean_kendall$P.Value, two = as.logical(result_wmean_kendall$two), invert = as.logical(result_wmean_kendall$invert_A))
        result_wmean_kendall$P.Value.I=two2one(result_wmean_kendall$P.Value, two = as.logical(result_wmean_kendall$two), invert = as.logical(result_wmean_kendall$invert_I))
        result_wmean_kendall_A=result_wmean_kendall[,c("geneset_id","P.Value.A")]
        result_wmean_kendall_I=result_wmean_kendall[,c("geneset_id","P.Value.I")]
        colnames(result_wmean_kendall_A)[2]=paste0(vector[i],".wmean_kendall")
        colnames(result_wmean_kendall_I)[2]=paste0(vector[i],".wmean_kendall")
        rm(result_wmean_kendall)
      }
      if("wmean_spearman" %in% enrichment_method){
        print("wmean_spearman will start.")
        result_wmean_spearman=spearman_test(geneset=geneset,value=value)
        result_wmean_spearman$two=TRUE
        result_wmean_spearman$invert_A=ifelse(result_wmean_spearman$r<=0,TRUE,FALSE)
        result_wmean_spearman$invert_I=ifelse(result_wmean_spearman$r>0,TRUE,FALSE)
        result_wmean_spearman$P.Value.A=two2one(result_wmean_spearman$P.Value, two = as.logical(result_wmean_spearman$two), invert = as.logical(result_wmean_spearman$invert_A))
        result_wmean_spearman$P.Value.I=two2one(result_wmean_spearman$P.Value, two = as.logical(result_wmean_spearman$two), invert = as.logical(result_wmean_spearman$invert_I))
        result_wmean_spearman_A=result_wmean_spearman[,c("geneset_id","P.Value.A")]
        result_wmean_spearman_I=result_wmean_spearman[,c("geneset_id","P.Value.I")]
        colnames(result_wmean_spearman_A)[2]=paste0(vector[i],".wmean_spearman")
        colnames(result_wmean_spearman_I)[2]=paste0(vector[i],".wmean_spearman")
        rm(result_wmean_spearman)
      }
      if("wmean_lm" %in% enrichment_method){
        print("wmean_lm will start.")
        result_wmean_lm=lm_test(geneset=geneset,value=value)
        result_wmean_lm$two=TRUE
        result_wmean_lm$invert_A=ifelse(result_wmean_lm$r<=0,TRUE,FALSE)
        result_wmean_lm$invert_I=ifelse(result_wmean_lm$r>0,TRUE,FALSE)
        result_wmean_lm$P.Value.A=two2one(result_wmean_lm$P.Value, two = as.logical(result_wmean_lm$two), invert = as.logical(result_wmean_lm$invert_A))
        result_wmean_lm$P.Value.I=two2one(result_wmean_lm$P.Value, two = as.logical(result_wmean_lm$two), invert = as.logical(result_wmean_lm$invert_I))
        result_wmean_lm_A=result_wmean_lm[,c("geneset_id","P.Value.A")]
        result_wmean_lm_I=result_wmean_lm[,c("geneset_id","P.Value.I")]
        colnames(result_wmean_lm_A)[2]=paste0(vector[i],".wmean_lm")
        colnames(result_wmean_lm_I)[2]=paste0(vector[i],".wmean_lm")
        rm(result_wmean_lm)
      }
      if("wmean_biweight" %in% enrichment_method){
        print("wmean_biweight will start.")
        result_wmean_biweight=biweight_test(geneset=geneset,value=value)
        result_wmean_biweight$two=TRUE
        result_wmean_biweight$invert_A=ifelse(result_wmean_biweight$r<=0,TRUE,FALSE)
        result_wmean_biweight$invert_I=ifelse(result_wmean_biweight$r>0,TRUE,FALSE)
        result_wmean_biweight$P.Value.A=two2one(result_wmean_biweight$P.Value, two = as.logical(result_wmean_biweight$two), invert = as.logical(result_wmean_biweight$invert_A))
        result_wmean_biweight$P.Value.I=two2one(result_wmean_biweight$P.Value, two = as.logical(result_wmean_biweight$two), invert = as.logical(result_wmean_biweight$invert_I))
        result_wmean_biweight_A=result_wmean_biweight[,c("geneset_id","P.Value.A")]
        result_wmean_biweight_I=result_wmean_biweight[,c("geneset_id","P.Value.I")]
        colnames(result_wmean_biweight_A)[2]=paste0(vector[i],".wmean_biweight")
        colnames(result_wmean_biweight_I)[2]=paste0(vector[i],".wmean_biweight")
        rm(result_wmean_biweight)
      }
      if("wmean_distance" %in% enrichment_method){
        print("wmean_distance will start.")
        result_wmean_distance=distance_test(geneset=geneset,value=value)
        result_wmean_distance$two=TRUE
        result_wmean_distance$invert_A=ifelse(result_wmean_distance$r<=0,TRUE,FALSE)
        result_wmean_distance$invert_I=ifelse(result_wmean_distance$r>0,TRUE,FALSE)
        result_wmean_distance$P.Value.A=two2one(result_wmean_distance$P.Value, two = as.logical(result_wmean_distance$two), invert = as.logical(result_wmean_distance$invert_A))
        result_wmean_distance$P.Value.I=two2one(result_wmean_distance$P.Value, two = as.logical(result_wmean_distance$two), invert = as.logical(result_wmean_distance$invert_I))
        result_wmean_distance_A=result_wmean_distance[,c("geneset_id","P.Value.A")]
        result_wmean_distance_I=result_wmean_distance[,c("geneset_id","P.Value.I")]
        colnames(result_wmean_distance_A)[2]=paste0(vector[i],".wmean_distance")
        colnames(result_wmean_distance_I)[2]=paste0(vector[i],".wmean_distance")
        rm(result_wmean_distance)
      }
      if("wmean_percentage" %in% enrichment_method){
        print("wmean_percentage will start.")
        result_wmean_percentage=percentage_test(geneset=geneset,value=value)
        result_wmean_percentage$two=TRUE
        result_wmean_percentage$invert_A=ifelse(result_wmean_percentage$r<=0,TRUE,FALSE)
        result_wmean_percentage$invert_I=ifelse(result_wmean_percentage$r>0,TRUE,FALSE)
        result_wmean_percentage$P.Value.A=two2one(result_wmean_percentage$P.Value, two = as.logical(result_wmean_percentage$two), invert = as.logical(result_wmean_percentage$invert_A))
        result_wmean_percentage$P.Value.I=two2one(result_wmean_percentage$P.Value, two = as.logical(result_wmean_percentage$two), invert = as.logical(result_wmean_percentage$invert_I))
        result_wmean_percentage_A=result_wmean_percentage[,c("geneset_id","P.Value.A")]
        result_wmean_percentage_I=result_wmean_percentage[,c("geneset_id","P.Value.I")]
        colnames(result_wmean_percentage_A)[2]=paste0(vector[i],".wmean_percentage")
        colnames(result_wmean_percentage_I)[2]=paste0(vector[i],".wmean_percentage")
        rm(result_wmean_percentage)
      }
      if("wmean_blomqvist" %in% enrichment_method){
        print("wmean_blomqvist will start.")
        result_wmean_blomqvist=blomqvist_test(geneset=geneset,value=value)
        result_wmean_blomqvist$two=TRUE
        result_wmean_blomqvist$invert_A=ifelse(result_wmean_blomqvist$r<=0,TRUE,FALSE)
        result_wmean_blomqvist$invert_I=ifelse(result_wmean_blomqvist$r>0,TRUE,FALSE)
        result_wmean_blomqvist$P.Value.A=two2one(result_wmean_blomqvist$P.Value, two = as.logical(result_wmean_blomqvist$two), invert = as.logical(result_wmean_blomqvist$invert_A))
        result_wmean_blomqvist$P.Value.I=two2one(result_wmean_blomqvist$P.Value, two = as.logical(result_wmean_blomqvist$two), invert = as.logical(result_wmean_blomqvist$invert_I))
        result_wmean_blomqvist_A=result_wmean_blomqvist[,c("geneset_id","P.Value.A")]
        result_wmean_blomqvist_I=result_wmean_blomqvist[,c("geneset_id","P.Value.I")]
        colnames(result_wmean_blomqvist_A)[2]=paste0(vector[i],".wmean_blomqvist")
        colnames(result_wmean_blomqvist_I)[2]=paste0(vector[i],".wmean_blomqvist")
        rm(result_wmean_blomqvist)
      }
      if("wmean_hoeffding" %in% enrichment_method){
        print("wmean_hoeffding will start.")
        result_wmean_hoeffding=hoeffding_test(geneset=geneset,value=value)
        result_wmean_hoeffding$two=TRUE
        result_wmean_hoeffding$invert_A=ifelse(result_wmean_hoeffding$r<=0,TRUE,FALSE)
        result_wmean_hoeffding$invert_I=ifelse(result_wmean_hoeffding$r>0,TRUE,FALSE)
        result_wmean_hoeffding$P.Value.A=two2one(result_wmean_hoeffding$P.Value, two = as.logical(result_wmean_hoeffding$two), invert = as.logical(result_wmean_hoeffding$invert_A))
        result_wmean_hoeffding$P.Value.I=two2one(result_wmean_hoeffding$P.Value, two = as.logical(result_wmean_hoeffding$two), invert = as.logical(result_wmean_hoeffding$invert_I))
        result_wmean_hoeffding_A=result_wmean_hoeffding[,c("geneset_id","P.Value.A")]
        result_wmean_hoeffding_I=result_wmean_hoeffding[,c("geneset_id","P.Value.I")]
        colnames(result_wmean_hoeffding_A)[2]=paste0(vector[i],".wmean_hoeffding")
        colnames(result_wmean_hoeffding_I)[2]=paste0(vector[i],".wmean_hoeffding")
        rm(result_wmean_hoeffding)
      }
      if("wmean_gamma" %in% enrichment_method){
        print("wmean_gamma will start.")
        result_wmean_gamma=gamma_test(geneset=geneset,value=value)
        result_wmean_gamma$two=TRUE
        result_wmean_gamma$invert_A=ifelse(result_wmean_gamma$r<=0,TRUE,FALSE)
        result_wmean_gamma$invert_I=ifelse(result_wmean_gamma$r>0,TRUE,FALSE)
        result_wmean_gamma$P.Value.A=two2one(result_wmean_gamma$P.Value, two = as.logical(result_wmean_gamma$two), invert = as.logical(result_wmean_gamma$invert_A))
        result_wmean_gamma$P.Value.I=two2one(result_wmean_gamma$P.Value, two = as.logical(result_wmean_gamma$two), invert = as.logical(result_wmean_gamma$invert_I))
        result_wmean_gamma_A=result_wmean_gamma[,c("geneset_id","P.Value.A")]
        result_wmean_gamma_I=result_wmean_gamma[,c("geneset_id","P.Value.I")]
        colnames(result_wmean_gamma_A)[2]=paste0(vector[i],".wmean_gamma")
        colnames(result_wmean_gamma_I)[2]=paste0(vector[i],".wmean_gamma")
        rm(result_wmean_gamma)
      }
    }else {
      print("Not using wmean.")
    }
    ##########################################norm_wmean
    methods_to_check=c("norm_wmean_pearson","norm_wmean_kendall","norm_wmean_spearman","norm_wmean_lm","norm_wmean_biweight","norm_wmean_distance","norm_wmean_percentage","norm_wmean_blomqvist","norm_wmean_hoeffding","norm_wmean_gamma")
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
      if("norm_wmean_pearson" %in% enrichment_method){
        print("norm_wmean_pearson will start.")
        result_norm_wmean_pearson=pearson_test(geneset=geneset,value=value)
        result_norm_wmean_pearson$two=TRUE
        result_norm_wmean_pearson$invert_A=ifelse(result_norm_wmean_pearson$r<=0,TRUE,FALSE)
        result_norm_wmean_pearson$invert_I=ifelse(result_norm_wmean_pearson$r>0,TRUE,FALSE)
        result_norm_wmean_pearson$P.Value.A=two2one(result_norm_wmean_pearson$P.Value, two = as.logical(result_norm_wmean_pearson$two), invert = as.logical(result_norm_wmean_pearson$invert_A))
        result_norm_wmean_pearson$P.Value.I=two2one(result_norm_wmean_pearson$P.Value, two = as.logical(result_norm_wmean_pearson$two), invert = as.logical(result_norm_wmean_pearson$invert_I))
        result_norm_wmean_pearson_A=result_norm_wmean_pearson[,c("geneset_id","P.Value.A")]
        result_norm_wmean_pearson_I=result_norm_wmean_pearson[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wmean_pearson_A)[2]=paste0(vector[i],".norm_wmean_pearson")
        colnames(result_norm_wmean_pearson_I)[2]=paste0(vector[i],".norm_wmean_pearson")
        rm(result_norm_wmean_pearson)
      }
      if("norm_wmean_kendall" %in% enrichment_method){
        print("norm_wmean_kendall will start.")
        result_norm_wmean_kendall=kendall_test(geneset=geneset,value=value)
        result_norm_wmean_kendall$two=TRUE
        result_norm_wmean_kendall$invert_A=ifelse(result_norm_wmean_kendall$r<=0,TRUE,FALSE)
        result_norm_wmean_kendall$invert_I=ifelse(result_norm_wmean_kendall$r>0,TRUE,FALSE)
        result_norm_wmean_kendall$P.Value.A=two2one(result_norm_wmean_kendall$P.Value, two = as.logical(result_norm_wmean_kendall$two), invert = as.logical(result_norm_wmean_kendall$invert_A))
        result_norm_wmean_kendall$P.Value.I=two2one(result_norm_wmean_kendall$P.Value, two = as.logical(result_norm_wmean_kendall$two), invert = as.logical(result_norm_wmean_kendall$invert_I))
        result_norm_wmean_kendall_A=result_norm_wmean_kendall[,c("geneset_id","P.Value.A")]
        result_norm_wmean_kendall_I=result_norm_wmean_kendall[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wmean_kendall_A)[2]=paste0(vector[i],".norm_wmean_kendall")
        colnames(result_norm_wmean_kendall_I)[2]=paste0(vector[i],".norm_wmean_kendall")
        rm(result_norm_wmean_kendall)
      }
      if("norm_wmean_spearman" %in% enrichment_method){
        print("norm_wmean_spearman will start.")
        result_norm_wmean_spearman=spearman_test(geneset=geneset,value=value)
        result_norm_wmean_spearman$two=TRUE
        result_norm_wmean_spearman$invert_A=ifelse(result_norm_wmean_spearman$r<=0,TRUE,FALSE)
        result_norm_wmean_spearman$invert_I=ifelse(result_norm_wmean_spearman$r>0,TRUE,FALSE)
        result_norm_wmean_spearman$P.Value.A=two2one(result_norm_wmean_spearman$P.Value, two = as.logical(result_norm_wmean_spearman$two), invert = as.logical(result_norm_wmean_spearman$invert_A))
        result_norm_wmean_spearman$P.Value.I=two2one(result_norm_wmean_spearman$P.Value, two = as.logical(result_norm_wmean_spearman$two), invert = as.logical(result_norm_wmean_spearman$invert_I))
        result_norm_wmean_spearman_A=result_norm_wmean_spearman[,c("geneset_id","P.Value.A")]
        result_norm_wmean_spearman_I=result_norm_wmean_spearman[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wmean_spearman_A)[2]=paste0(vector[i],".norm_wmean_spearman")
        colnames(result_norm_wmean_spearman_I)[2]=paste0(vector[i],".norm_wmean_spearman")
        rm(result_norm_wmean_spearman)
      }
      if("norm_wmean_lm" %in% enrichment_method){
        print("norm_wmean_lm will start.")
        result_norm_wmean_lm=lm_test(geneset=geneset,value=value)
        result_norm_wmean_lm$two=TRUE
        result_norm_wmean_lm$invert_A=ifelse(result_norm_wmean_lm$r<=0,TRUE,FALSE)
        result_norm_wmean_lm$invert_I=ifelse(result_norm_wmean_lm$r>0,TRUE,FALSE)
        result_norm_wmean_lm$P.Value.A=two2one(result_norm_wmean_lm$P.Value, two = as.logical(result_norm_wmean_lm$two), invert = as.logical(result_norm_wmean_lm$invert_A))
        result_norm_wmean_lm$P.Value.I=two2one(result_norm_wmean_lm$P.Value, two = as.logical(result_norm_wmean_lm$two), invert = as.logical(result_norm_wmean_lm$invert_I))
        result_norm_wmean_lm_A=result_norm_wmean_lm[,c("geneset_id","P.Value.A")]
        result_norm_wmean_lm_I=result_norm_wmean_lm[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wmean_lm_A)[2]=paste0(vector[i],".norm_wmean_lm")
        colnames(result_norm_wmean_lm_I)[2]=paste0(vector[i],".norm_wmean_lm")
        rm(result_norm_wmean_lm)
      }
      if("norm_wmean_biweight" %in% enrichment_method){
        print("norm_wmean_biweight will start.")
        result_norm_wmean_biweight=biweight_test(geneset=geneset,value=value)
        result_norm_wmean_biweight$two=TRUE
        result_norm_wmean_biweight$invert_A=ifelse(result_norm_wmean_biweight$r<=0,TRUE,FALSE)
        result_norm_wmean_biweight$invert_I=ifelse(result_norm_wmean_biweight$r>0,TRUE,FALSE)
        result_norm_wmean_biweight$P.Value.A=two2one(result_norm_wmean_biweight$P.Value, two = as.logical(result_norm_wmean_biweight$two), invert = as.logical(result_norm_wmean_biweight$invert_A))
        result_norm_wmean_biweight$P.Value.I=two2one(result_norm_wmean_biweight$P.Value, two = as.logical(result_norm_wmean_biweight$two), invert = as.logical(result_norm_wmean_biweight$invert_I))
        result_norm_wmean_biweight_A=result_norm_wmean_biweight[,c("geneset_id","P.Value.A")]
        result_norm_wmean_biweight_I=result_norm_wmean_biweight[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wmean_biweight_A)[2]=paste0(vector[i],".norm_wmean_biweight")
        colnames(result_norm_wmean_biweight_I)[2]=paste0(vector[i],".norm_wmean_biweight")
        rm(result_norm_wmean_biweight)
      }
      if("norm_wmean_distance" %in% enrichment_method){
        print("norm_wmean_distance will start.")
        result_norm_wmean_distance=distance_test(geneset=geneset,value=value)
        result_norm_wmean_distance$two=TRUE
        result_norm_wmean_distance$invert_A=ifelse(result_norm_wmean_distance$r<=0,TRUE,FALSE)
        result_norm_wmean_distance$invert_I=ifelse(result_norm_wmean_distance$r>0,TRUE,FALSE)
        result_norm_wmean_distance$P.Value.A=two2one(result_norm_wmean_distance$P.Value, two = as.logical(result_norm_wmean_distance$two), invert = as.logical(result_norm_wmean_distance$invert_A))
        result_norm_wmean_distance$P.Value.I=two2one(result_norm_wmean_distance$P.Value, two = as.logical(result_norm_wmean_distance$two), invert = as.logical(result_norm_wmean_distance$invert_I))
        result_norm_wmean_distance_A=result_norm_wmean_distance[,c("geneset_id","P.Value.A")]
        result_norm_wmean_distance_I=result_norm_wmean_distance[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wmean_distance_A)[2]=paste0(vector[i],".norm_wmean_distance")
        colnames(result_norm_wmean_distance_I)[2]=paste0(vector[i],".norm_wmean_distance")
        rm(result_norm_wmean_distance)
      }
      if("norm_wmean_percentage" %in% enrichment_method){
        print("norm_wmean_percentage will start.")
        result_norm_wmean_percentage=percentage_test(geneset=geneset,value=value)
        result_norm_wmean_percentage$two=TRUE
        result_norm_wmean_percentage$invert_A=ifelse(result_norm_wmean_percentage$r<=0,TRUE,FALSE)
        result_norm_wmean_percentage$invert_I=ifelse(result_norm_wmean_percentage$r>0,TRUE,FALSE)
        result_norm_wmean_percentage$P.Value.A=two2one(result_norm_wmean_percentage$P.Value, two = as.logical(result_norm_wmean_percentage$two), invert = as.logical(result_norm_wmean_percentage$invert_A))
        result_norm_wmean_percentage$P.Value.I=two2one(result_norm_wmean_percentage$P.Value, two = as.logical(result_norm_wmean_percentage$two), invert = as.logical(result_norm_wmean_percentage$invert_I))
        result_norm_wmean_percentage_A=result_norm_wmean_percentage[,c("geneset_id","P.Value.A")]
        result_norm_wmean_percentage_I=result_norm_wmean_percentage[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wmean_percentage_A)[2]=paste0(vector[i],".norm_wmean_percentage")
        colnames(result_norm_wmean_percentage_I)[2]=paste0(vector[i],".norm_wmean_percentage")
        rm(result_norm_wmean_percentage)
      }
      if("norm_wmean_blomqvist" %in% enrichment_method){
        print("norm_wmean_blomqvist will start.")
        result_norm_wmean_blomqvist=blomqvist_test(geneset=geneset,value=value)
        result_norm_wmean_blomqvist$two=TRUE
        result_norm_wmean_blomqvist$invert_A=ifelse(result_norm_wmean_blomqvist$r<=0,TRUE,FALSE)
        result_norm_wmean_blomqvist$invert_I=ifelse(result_norm_wmean_blomqvist$r>0,TRUE,FALSE)
        result_norm_wmean_blomqvist$P.Value.A=two2one(result_norm_wmean_blomqvist$P.Value, two = as.logical(result_norm_wmean_blomqvist$two), invert = as.logical(result_norm_wmean_blomqvist$invert_A))
        result_norm_wmean_blomqvist$P.Value.I=two2one(result_norm_wmean_blomqvist$P.Value, two = as.logical(result_norm_wmean_blomqvist$two), invert = as.logical(result_norm_wmean_blomqvist$invert_I))
        result_norm_wmean_blomqvist_A=result_norm_wmean_blomqvist[,c("geneset_id","P.Value.A")]
        result_norm_wmean_blomqvist_I=result_norm_wmean_blomqvist[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wmean_blomqvist_A)[2]=paste0(vector[i],".norm_wmean_blomqvist")
        colnames(result_norm_wmean_blomqvist_I)[2]=paste0(vector[i],".norm_wmean_blomqvist")
        rm(result_norm_wmean_blomqvist)
      }
      if("norm_wmean_hoeffding" %in% enrichment_method){
        print("norm_wmean_hoeffding will start.")
        result_norm_wmean_hoeffding=hoeffding_test(geneset=geneset,value=value)
        result_norm_wmean_hoeffding$two=TRUE
        result_norm_wmean_hoeffding$invert_A=ifelse(result_norm_wmean_hoeffding$r<=0,TRUE,FALSE)
        result_norm_wmean_hoeffding$invert_I=ifelse(result_norm_wmean_hoeffding$r>0,TRUE,FALSE)
        result_norm_wmean_hoeffding$P.Value.A=two2one(result_norm_wmean_hoeffding$P.Value, two = as.logical(result_norm_wmean_hoeffding$two), invert = as.logical(result_norm_wmean_hoeffding$invert_A))
        result_norm_wmean_hoeffding$P.Value.I=two2one(result_norm_wmean_hoeffding$P.Value, two = as.logical(result_norm_wmean_hoeffding$two), invert = as.logical(result_norm_wmean_hoeffding$invert_I))
        result_norm_wmean_hoeffding_A=result_norm_wmean_hoeffding[,c("geneset_id","P.Value.A")]
        result_norm_wmean_hoeffding_I=result_norm_wmean_hoeffding[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wmean_hoeffding_A)[2]=paste0(vector[i],".norm_wmean_hoeffding")
        colnames(result_norm_wmean_hoeffding_I)[2]=paste0(vector[i],".norm_wmean_hoeffding")
        rm(result_norm_wmean_hoeffding)
      }
      if("norm_wmean_gamma" %in% enrichment_method){
        print("norm_wmean_gamma will start.")
        result_norm_wmean_gamma=gamma_test(geneset=geneset,value=value)
        result_norm_wmean_gamma$two=TRUE
        result_norm_wmean_gamma$invert_A=ifelse(result_norm_wmean_gamma$r<=0,TRUE,FALSE)
        result_norm_wmean_gamma$invert_I=ifelse(result_norm_wmean_gamma$r>0,TRUE,FALSE)
        result_norm_wmean_gamma$P.Value.A=two2one(result_norm_wmean_gamma$P.Value, two = as.logical(result_norm_wmean_gamma$two), invert = as.logical(result_norm_wmean_gamma$invert_A))
        result_norm_wmean_gamma$P.Value.I=two2one(result_norm_wmean_gamma$P.Value, two = as.logical(result_norm_wmean_gamma$two), invert = as.logical(result_norm_wmean_gamma$invert_I))
        result_norm_wmean_gamma_A=result_norm_wmean_gamma[,c("geneset_id","P.Value.A")]
        result_norm_wmean_gamma_I=result_norm_wmean_gamma[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wmean_gamma_A)[2]=paste0(vector[i],".norm_wmean_gamma")
        colnames(result_norm_wmean_gamma_I)[2]=paste0(vector[i],".norm_wmean_gamma")
        rm(result_norm_wmean_gamma)
      }
    }else {
      print("Not using norm_wmean.")
    }
    ##########################################corr_wmean
    methods_to_check=c("corr_wmean_pearson","corr_wmean_kendall","corr_wmean_spearman","corr_wmean_lm","corr_wmean_biweight","corr_wmean_distance","corr_wmean_percentage","corr_wmean_blomqvist","corr_wmean_hoeffding","corr_wmean_gamma")
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
      if("corr_wmean_pearson" %in% enrichment_method){
        print("corr_wmean_pearson will start.")
        result_corr_wmean_pearson=pearson_test(geneset=geneset,value=value)
        result_corr_wmean_pearson$two=TRUE
        result_corr_wmean_pearson$invert_A=ifelse(result_corr_wmean_pearson$r<=0,TRUE,FALSE)
        result_corr_wmean_pearson$invert_I=ifelse(result_corr_wmean_pearson$r>0,TRUE,FALSE)
        result_corr_wmean_pearson$P.Value.A=two2one(result_corr_wmean_pearson$P.Value, two = as.logical(result_corr_wmean_pearson$two), invert = as.logical(result_corr_wmean_pearson$invert_A))
        result_corr_wmean_pearson$P.Value.I=two2one(result_corr_wmean_pearson$P.Value, two = as.logical(result_corr_wmean_pearson$two), invert = as.logical(result_corr_wmean_pearson$invert_I))
        result_corr_wmean_pearson_A=result_corr_wmean_pearson[,c("geneset_id","P.Value.A")]
        result_corr_wmean_pearson_I=result_corr_wmean_pearson[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wmean_pearson_A)[2]=paste0(vector[i],".corr_wmean_pearson")
        colnames(result_corr_wmean_pearson_I)[2]=paste0(vector[i],".corr_wmean_pearson")
        rm(result_corr_wmean_pearson)
      }
      if("corr_wmean_kendall" %in% enrichment_method){
        print("corr_wmean_kendall will start.")
        result_corr_wmean_kendall=kendall_test(geneset=geneset,value=value)
        result_corr_wmean_kendall$two=TRUE
        result_corr_wmean_kendall$invert_A=ifelse(result_corr_wmean_kendall$r<=0,TRUE,FALSE)
        result_corr_wmean_kendall$invert_I=ifelse(result_corr_wmean_kendall$r>0,TRUE,FALSE)
        result_corr_wmean_kendall$P.Value.A=two2one(result_corr_wmean_kendall$P.Value, two = as.logical(result_corr_wmean_kendall$two), invert = as.logical(result_corr_wmean_kendall$invert_A))
        result_corr_wmean_kendall$P.Value.I=two2one(result_corr_wmean_kendall$P.Value, two = as.logical(result_corr_wmean_kendall$two), invert = as.logical(result_corr_wmean_kendall$invert_I))
        result_corr_wmean_kendall_A=result_corr_wmean_kendall[,c("geneset_id","P.Value.A")]
        result_corr_wmean_kendall_I=result_corr_wmean_kendall[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wmean_kendall_A)[2]=paste0(vector[i],".corr_wmean_kendall")
        colnames(result_corr_wmean_kendall_I)[2]=paste0(vector[i],".corr_wmean_kendall")
        rm(result_corr_wmean_kendall)
      }
      if("corr_wmean_spearman" %in% enrichment_method){
        print("corr_wmean_spearman will start.")
        result_corr_wmean_spearman=spearman_test(geneset=geneset,value=value)
        result_corr_wmean_spearman$two=TRUE
        result_corr_wmean_spearman$invert_A=ifelse(result_corr_wmean_spearman$r<=0,TRUE,FALSE)
        result_corr_wmean_spearman$invert_I=ifelse(result_corr_wmean_spearman$r>0,TRUE,FALSE)
        result_corr_wmean_spearman$P.Value.A=two2one(result_corr_wmean_spearman$P.Value, two = as.logical(result_corr_wmean_spearman$two), invert = as.logical(result_corr_wmean_spearman$invert_A))
        result_corr_wmean_spearman$P.Value.I=two2one(result_corr_wmean_spearman$P.Value, two = as.logical(result_corr_wmean_spearman$two), invert = as.logical(result_corr_wmean_spearman$invert_I))
        result_corr_wmean_spearman_A=result_corr_wmean_spearman[,c("geneset_id","P.Value.A")]
        result_corr_wmean_spearman_I=result_corr_wmean_spearman[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wmean_spearman_A)[2]=paste0(vector[i],".corr_wmean_spearman")
        colnames(result_corr_wmean_spearman_I)[2]=paste0(vector[i],".corr_wmean_spearman")
        rm(result_corr_wmean_spearman)
      }
      if("corr_wmean_lm" %in% enrichment_method){
        print("corr_wmean_lm will start.")
        result_corr_wmean_lm=lm_test(geneset=geneset,value=value)
        result_corr_wmean_lm$two=TRUE
        result_corr_wmean_lm$invert_A=ifelse(result_corr_wmean_lm$r<=0,TRUE,FALSE)
        result_corr_wmean_lm$invert_I=ifelse(result_corr_wmean_lm$r>0,TRUE,FALSE)
        result_corr_wmean_lm$P.Value.A=two2one(result_corr_wmean_lm$P.Value, two = as.logical(result_corr_wmean_lm$two), invert = as.logical(result_corr_wmean_lm$invert_A))
        result_corr_wmean_lm$P.Value.I=two2one(result_corr_wmean_lm$P.Value, two = as.logical(result_corr_wmean_lm$two), invert = as.logical(result_corr_wmean_lm$invert_I))
        result_corr_wmean_lm_A=result_corr_wmean_lm[,c("geneset_id","P.Value.A")]
        result_corr_wmean_lm_I=result_corr_wmean_lm[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wmean_lm_A)[2]=paste0(vector[i],".corr_wmean_lm")
        colnames(result_corr_wmean_lm_I)[2]=paste0(vector[i],".corr_wmean_lm")
        rm(result_corr_wmean_lm)
      }
      if("corr_wmean_biweight" %in% enrichment_method){
        print("corr_wmean_biweight will start.")
        result_corr_wmean_biweight=biweight_test(geneset=geneset,value=value)
        result_corr_wmean_biweight$two=TRUE
        result_corr_wmean_biweight$invert_A=ifelse(result_corr_wmean_biweight$r<=0,TRUE,FALSE)
        result_corr_wmean_biweight$invert_I=ifelse(result_corr_wmean_biweight$r>0,TRUE,FALSE)
        result_corr_wmean_biweight$P.Value.A=two2one(result_corr_wmean_biweight$P.Value, two = as.logical(result_corr_wmean_biweight$two), invert = as.logical(result_corr_wmean_biweight$invert_A))
        result_corr_wmean_biweight$P.Value.I=two2one(result_corr_wmean_biweight$P.Value, two = as.logical(result_corr_wmean_biweight$two), invert = as.logical(result_corr_wmean_biweight$invert_I))
        result_corr_wmean_biweight_A=result_corr_wmean_biweight[,c("geneset_id","P.Value.A")]
        result_corr_wmean_biweight_I=result_corr_wmean_biweight[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wmean_biweight_A)[2]=paste0(vector[i],".corr_wmean_biweight")
        colnames(result_corr_wmean_biweight_I)[2]=paste0(vector[i],".corr_wmean_biweight")
        rm(result_corr_wmean_biweight)
      }
      if("corr_wmean_distance" %in% enrichment_method){
        print("corr_wmean_distance will start.")
        result_corr_wmean_distance=distance_test(geneset=geneset,value=value)
        result_corr_wmean_distance$two=TRUE
        result_corr_wmean_distance$invert_A=ifelse(result_corr_wmean_distance$r<=0,TRUE,FALSE)
        result_corr_wmean_distance$invert_I=ifelse(result_corr_wmean_distance$r>0,TRUE,FALSE)
        result_corr_wmean_distance$P.Value.A=two2one(result_corr_wmean_distance$P.Value, two = as.logical(result_corr_wmean_distance$two), invert = as.logical(result_corr_wmean_distance$invert_A))
        result_corr_wmean_distance$P.Value.I=two2one(result_corr_wmean_distance$P.Value, two = as.logical(result_corr_wmean_distance$two), invert = as.logical(result_corr_wmean_distance$invert_I))
        result_corr_wmean_distance_A=result_corr_wmean_distance[,c("geneset_id","P.Value.A")]
        result_corr_wmean_distance_I=result_corr_wmean_distance[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wmean_distance_A)[2]=paste0(vector[i],".corr_wmean_distance")
        colnames(result_corr_wmean_distance_I)[2]=paste0(vector[i],".corr_wmean_distance")
        rm(result_corr_wmean_distance)
      }
      if("corr_wmean_percentage" %in% enrichment_method){
        print("corr_wmean_percentage will start.")
        result_corr_wmean_percentage=percentage_test(geneset=geneset,value=value)
        result_corr_wmean_percentage$two=TRUE
        result_corr_wmean_percentage$invert_A=ifelse(result_corr_wmean_percentage$r<=0,TRUE,FALSE)
        result_corr_wmean_percentage$invert_I=ifelse(result_corr_wmean_percentage$r>0,TRUE,FALSE)
        result_corr_wmean_percentage$P.Value.A=two2one(result_corr_wmean_percentage$P.Value, two = as.logical(result_corr_wmean_percentage$two), invert = as.logical(result_corr_wmean_percentage$invert_A))
        result_corr_wmean_percentage$P.Value.I=two2one(result_corr_wmean_percentage$P.Value, two = as.logical(result_corr_wmean_percentage$two), invert = as.logical(result_corr_wmean_percentage$invert_I))
        result_corr_wmean_percentage_A=result_corr_wmean_percentage[,c("geneset_id","P.Value.A")]
        result_corr_wmean_percentage_I=result_corr_wmean_percentage[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wmean_percentage_A)[2]=paste0(vector[i],".corr_wmean_percentage")
        colnames(result_corr_wmean_percentage_I)[2]=paste0(vector[i],".corr_wmean_percentage")
        rm(result_corr_wmean_percentage)
      }
      if("corr_wmean_blomqvist" %in% enrichment_method){
        print("corr_wmean_blomqvist will start.")
        result_corr_wmean_blomqvist=blomqvist_test(geneset=geneset,value=value)
        result_corr_wmean_blomqvist$two=TRUE
        result_corr_wmean_blomqvist$invert_A=ifelse(result_corr_wmean_blomqvist$r<=0,TRUE,FALSE)
        result_corr_wmean_blomqvist$invert_I=ifelse(result_corr_wmean_blomqvist$r>0,TRUE,FALSE)
        result_corr_wmean_blomqvist$P.Value.A=two2one(result_corr_wmean_blomqvist$P.Value, two = as.logical(result_corr_wmean_blomqvist$two), invert = as.logical(result_corr_wmean_blomqvist$invert_A))
        result_corr_wmean_blomqvist$P.Value.I=two2one(result_corr_wmean_blomqvist$P.Value, two = as.logical(result_corr_wmean_blomqvist$two), invert = as.logical(result_corr_wmean_blomqvist$invert_I))
        result_corr_wmean_blomqvist_A=result_corr_wmean_blomqvist[,c("geneset_id","P.Value.A")]
        result_corr_wmean_blomqvist_I=result_corr_wmean_blomqvist[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wmean_blomqvist_A)[2]=paste0(vector[i],".corr_wmean_blomqvist")
        colnames(result_corr_wmean_blomqvist_I)[2]=paste0(vector[i],".corr_wmean_blomqvist")
        rm(result_corr_wmean_blomqvist)
      }
      if("corr_wmean_hoeffding" %in% enrichment_method){
        print("corr_wmean_hoeffding will start.")
        result_corr_wmean_hoeffding=hoeffding_test(geneset=geneset,value=value)
        result_corr_wmean_hoeffding$two=TRUE
        result_corr_wmean_hoeffding$invert_A=ifelse(result_corr_wmean_hoeffding$r<=0,TRUE,FALSE)
        result_corr_wmean_hoeffding$invert_I=ifelse(result_corr_wmean_hoeffding$r>0,TRUE,FALSE)
        result_corr_wmean_hoeffding$P.Value.A=two2one(result_corr_wmean_hoeffding$P.Value, two = as.logical(result_corr_wmean_hoeffding$two), invert = as.logical(result_corr_wmean_hoeffding$invert_A))
        result_corr_wmean_hoeffding$P.Value.I=two2one(result_corr_wmean_hoeffding$P.Value, two = as.logical(result_corr_wmean_hoeffding$two), invert = as.logical(result_corr_wmean_hoeffding$invert_I))
        result_corr_wmean_hoeffding_A=result_corr_wmean_hoeffding[,c("geneset_id","P.Value.A")]
        result_corr_wmean_hoeffding_I=result_corr_wmean_hoeffding[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wmean_hoeffding_A)[2]=paste0(vector[i],".corr_wmean_hoeffding")
        colnames(result_corr_wmean_hoeffding_I)[2]=paste0(vector[i],".corr_wmean_hoeffding")
        rm(result_corr_wmean_hoeffding)
      }
      if("corr_wmean_gamma" %in% enrichment_method){
        print("corr_wmean_gamma will start.")
        result_corr_wmean_gamma=gamma_test(geneset=geneset,value=value)
        result_corr_wmean_gamma$two=TRUE
        result_corr_wmean_gamma$invert_A=ifelse(result_corr_wmean_gamma$r<=0,TRUE,FALSE)
        result_corr_wmean_gamma$invert_I=ifelse(result_corr_wmean_gamma$r>0,TRUE,FALSE)
        result_corr_wmean_gamma$P.Value.A=two2one(result_corr_wmean_gamma$P.Value, two = as.logical(result_corr_wmean_gamma$two), invert = as.logical(result_corr_wmean_gamma$invert_A))
        result_corr_wmean_gamma$P.Value.I=two2one(result_corr_wmean_gamma$P.Value, two = as.logical(result_corr_wmean_gamma$two), invert = as.logical(result_corr_wmean_gamma$invert_I))
        result_corr_wmean_gamma_A=result_corr_wmean_gamma[,c("geneset_id","P.Value.A")]
        result_corr_wmean_gamma_I=result_corr_wmean_gamma[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wmean_gamma_A)[2]=paste0(vector[i],".corr_wmean_gamma")
        colnames(result_corr_wmean_gamma_I)[2]=paste0(vector[i],".corr_wmean_gamma")
        rm(result_corr_wmean_gamma)
      }
    }else {
      print("Not using corr_wmean.")
    }
    ##########################################wsum
    methods_to_check=c("wsum_pearson","wsum_kendall","wsum_spearman","wsum_lm","wsum_biweight","wsum_distance","wsum_percentage","wsum_blomqvist","wsum_hoeffding","wsum_gamma")
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
      if("wsum_pearson" %in% enrichment_method){
        print("wsum_pearson will start.")
        result_wsum_pearson=pearson_test(geneset=geneset,value=value)
        result_wsum_pearson$two=TRUE
        result_wsum_pearson$invert_A=ifelse(result_wsum_pearson$r<=0,TRUE,FALSE)
        result_wsum_pearson$invert_I=ifelse(result_wsum_pearson$r>0,TRUE,FALSE)
        result_wsum_pearson$P.Value.A=two2one(result_wsum_pearson$P.Value, two = as.logical(result_wsum_pearson$two), invert = as.logical(result_wsum_pearson$invert_A))
        result_wsum_pearson$P.Value.I=two2one(result_wsum_pearson$P.Value, two = as.logical(result_wsum_pearson$two), invert = as.logical(result_wsum_pearson$invert_I))
        result_wsum_pearson_A=result_wsum_pearson[,c("geneset_id","P.Value.A")]
        result_wsum_pearson_I=result_wsum_pearson[,c("geneset_id","P.Value.I")]
        colnames(result_wsum_pearson_A)[2]=paste0(vector[i],".wsum_pearson")
        colnames(result_wsum_pearson_I)[2]=paste0(vector[i],".wsum_pearson")
        rm(result_wsum_pearson)
      }
      if("wsum_kendall" %in% enrichment_method){
        print("wsum_kendall will start.")
        result_wsum_kendall=kendall_test(geneset=geneset,value=value)
        result_wsum_kendall$two=TRUE
        result_wsum_kendall$invert_A=ifelse(result_wsum_kendall$r<=0,TRUE,FALSE)
        result_wsum_kendall$invert_I=ifelse(result_wsum_kendall$r>0,TRUE,FALSE)
        result_wsum_kendall$P.Value.A=two2one(result_wsum_kendall$P.Value, two = as.logical(result_wsum_kendall$two), invert = as.logical(result_wsum_kendall$invert_A))
        result_wsum_kendall$P.Value.I=two2one(result_wsum_kendall$P.Value, two = as.logical(result_wsum_kendall$two), invert = as.logical(result_wsum_kendall$invert_I))
        result_wsum_kendall_A=result_wsum_kendall[,c("geneset_id","P.Value.A")]
        result_wsum_kendall_I=result_wsum_kendall[,c("geneset_id","P.Value.I")]
        colnames(result_wsum_kendall_A)[2]=paste0(vector[i],".wsum_kendall")
        colnames(result_wsum_kendall_I)[2]=paste0(vector[i],".wsum_kendall")
        rm(result_wsum_kendall)
      }
      if("wsum_spearman" %in% enrichment_method){
        print("wsum_spearman will start.")
        result_wsum_spearman=spearman_test(geneset=geneset,value=value)
        result_wsum_spearman$two=TRUE
        result_wsum_spearman$invert_A=ifelse(result_wsum_spearman$r<=0,TRUE,FALSE)
        result_wsum_spearman$invert_I=ifelse(result_wsum_spearman$r>0,TRUE,FALSE)
        result_wsum_spearman$P.Value.A=two2one(result_wsum_spearman$P.Value, two = as.logical(result_wsum_spearman$two), invert = as.logical(result_wsum_spearman$invert_A))
        result_wsum_spearman$P.Value.I=two2one(result_wsum_spearman$P.Value, two = as.logical(result_wsum_spearman$two), invert = as.logical(result_wsum_spearman$invert_I))
        result_wsum_spearman_A=result_wsum_spearman[,c("geneset_id","P.Value.A")]
        result_wsum_spearman_I=result_wsum_spearman[,c("geneset_id","P.Value.I")]
        colnames(result_wsum_spearman_A)[2]=paste0(vector[i],".wsum_spearman")
        colnames(result_wsum_spearman_I)[2]=paste0(vector[i],".wsum_spearman")
        rm(result_wsum_spearman)
      }
      if("wsum_lm" %in% enrichment_method){
        print("wsum_lm will start.")
        result_wsum_lm=lm_test(geneset=geneset,value=value)
        result_wsum_lm$two=TRUE
        result_wsum_lm$invert_A=ifelse(result_wsum_lm$r<=0,TRUE,FALSE)
        result_wsum_lm$invert_I=ifelse(result_wsum_lm$r>0,TRUE,FALSE)
        result_wsum_lm$P.Value.A=two2one(result_wsum_lm$P.Value, two = as.logical(result_wsum_lm$two), invert = as.logical(result_wsum_lm$invert_A))
        result_wsum_lm$P.Value.I=two2one(result_wsum_lm$P.Value, two = as.logical(result_wsum_lm$two), invert = as.logical(result_wsum_lm$invert_I))
        result_wsum_lm_A=result_wsum_lm[,c("geneset_id","P.Value.A")]
        result_wsum_lm_I=result_wsum_lm[,c("geneset_id","P.Value.I")]
        colnames(result_wsum_lm_A)[2]=paste0(vector[i],".wsum_lm")
        colnames(result_wsum_lm_I)[2]=paste0(vector[i],".wsum_lm")
        rm(result_wsum_lm)
      }
      if("wsum_biweight" %in% enrichment_method){
        print("wsum_biweight will start.")
        result_wsum_biweight=biweight_test(geneset=geneset,value=value)
        result_wsum_biweight$two=TRUE
        result_wsum_biweight$invert_A=ifelse(result_wsum_biweight$r<=0,TRUE,FALSE)
        result_wsum_biweight$invert_I=ifelse(result_wsum_biweight$r>0,TRUE,FALSE)
        result_wsum_biweight$P.Value.A=two2one(result_wsum_biweight$P.Value, two = as.logical(result_wsum_biweight$two), invert = as.logical(result_wsum_biweight$invert_A))
        result_wsum_biweight$P.Value.I=two2one(result_wsum_biweight$P.Value, two = as.logical(result_wsum_biweight$two), invert = as.logical(result_wsum_biweight$invert_I))
        result_wsum_biweight_A=result_wsum_biweight[,c("geneset_id","P.Value.A")]
        result_wsum_biweight_I=result_wsum_biweight[,c("geneset_id","P.Value.I")]
        colnames(result_wsum_biweight_A)[2]=paste0(vector[i],".wsum_biweight")
        colnames(result_wsum_biweight_I)[2]=paste0(vector[i],".wsum_biweight")
        rm(result_wsum_biweight)
      }
      if("wsum_distance" %in% enrichment_method){
        print("wsum_distance will start.")
        result_wsum_distance=distance_test(geneset=geneset,value=value)
        result_wsum_distance$two=TRUE
        result_wsum_distance$invert_A=ifelse(result_wsum_distance$r<=0,TRUE,FALSE)
        result_wsum_distance$invert_I=ifelse(result_wsum_distance$r>0,TRUE,FALSE)
        result_wsum_distance$P.Value.A=two2one(result_wsum_distance$P.Value, two = as.logical(result_wsum_distance$two), invert = as.logical(result_wsum_distance$invert_A))
        result_wsum_distance$P.Value.I=two2one(result_wsum_distance$P.Value, two = as.logical(result_wsum_distance$two), invert = as.logical(result_wsum_distance$invert_I))
        result_wsum_distance_A=result_wsum_distance[,c("geneset_id","P.Value.A")]
        result_wsum_distance_I=result_wsum_distance[,c("geneset_id","P.Value.I")]
        colnames(result_wsum_distance_A)[2]=paste0(vector[i],".wsum_distance")
        colnames(result_wsum_distance_I)[2]=paste0(vector[i],".wsum_distance")
        rm(result_wsum_distance)
      }
      if("wsum_percentage" %in% enrichment_method){
        print("wsum_percentage will start.")
        result_wsum_percentage=percentage_test(geneset=geneset,value=value)
        result_wsum_percentage$two=TRUE
        result_wsum_percentage$invert_A=ifelse(result_wsum_percentage$r<=0,TRUE,FALSE)
        result_wsum_percentage$invert_I=ifelse(result_wsum_percentage$r>0,TRUE,FALSE)
        result_wsum_percentage$P.Value.A=two2one(result_wsum_percentage$P.Value, two = as.logical(result_wsum_percentage$two), invert = as.logical(result_wsum_percentage$invert_A))
        result_wsum_percentage$P.Value.I=two2one(result_wsum_percentage$P.Value, two = as.logical(result_wsum_percentage$two), invert = as.logical(result_wsum_percentage$invert_I))
        result_wsum_percentage_A=result_wsum_percentage[,c("geneset_id","P.Value.A")]
        result_wsum_percentage_I=result_wsum_percentage[,c("geneset_id","P.Value.I")]
        colnames(result_wsum_percentage_A)[2]=paste0(vector[i],".wsum_percentage")
        colnames(result_wsum_percentage_I)[2]=paste0(vector[i],".wsum_percentage")
        rm(result_wsum_percentage)
      }
      if("wsum_blomqvist" %in% enrichment_method){
        print("wsum_blomqvist will start.")
        result_wsum_blomqvist=blomqvist_test(geneset=geneset,value=value)
        result_wsum_blomqvist$two=TRUE
        result_wsum_blomqvist$invert_A=ifelse(result_wsum_blomqvist$r<=0,TRUE,FALSE)
        result_wsum_blomqvist$invert_I=ifelse(result_wsum_blomqvist$r>0,TRUE,FALSE)
        result_wsum_blomqvist$P.Value.A=two2one(result_wsum_blomqvist$P.Value, two = as.logical(result_wsum_blomqvist$two), invert = as.logical(result_wsum_blomqvist$invert_A))
        result_wsum_blomqvist$P.Value.I=two2one(result_wsum_blomqvist$P.Value, two = as.logical(result_wsum_blomqvist$two), invert = as.logical(result_wsum_blomqvist$invert_I))
        result_wsum_blomqvist_A=result_wsum_blomqvist[,c("geneset_id","P.Value.A")]
        result_wsum_blomqvist_I=result_wsum_blomqvist[,c("geneset_id","P.Value.I")]
        colnames(result_wsum_blomqvist_A)[2]=paste0(vector[i],".wsum_blomqvist")
        colnames(result_wsum_blomqvist_I)[2]=paste0(vector[i],".wsum_blomqvist")
        rm(result_wsum_blomqvist)
      }
      if("wsum_hoeffding" %in% enrichment_method){
        print("wsum_hoeffding will start.")
        result_wsum_hoeffding=hoeffding_test(geneset=geneset,value=value)
        result_wsum_hoeffding$two=TRUE
        result_wsum_hoeffding$invert_A=ifelse(result_wsum_hoeffding$r<=0,TRUE,FALSE)
        result_wsum_hoeffding$invert_I=ifelse(result_wsum_hoeffding$r>0,TRUE,FALSE)
        result_wsum_hoeffding$P.Value.A=two2one(result_wsum_hoeffding$P.Value, two = as.logical(result_wsum_hoeffding$two), invert = as.logical(result_wsum_hoeffding$invert_A))
        result_wsum_hoeffding$P.Value.I=two2one(result_wsum_hoeffding$P.Value, two = as.logical(result_wsum_hoeffding$two), invert = as.logical(result_wsum_hoeffding$invert_I))
        result_wsum_hoeffding_A=result_wsum_hoeffding[,c("geneset_id","P.Value.A")]
        result_wsum_hoeffding_I=result_wsum_hoeffding[,c("geneset_id","P.Value.I")]
        colnames(result_wsum_hoeffding_A)[2]=paste0(vector[i],".wsum_hoeffding")
        colnames(result_wsum_hoeffding_I)[2]=paste0(vector[i],".wsum_hoeffding")
        rm(result_wsum_hoeffding)
      }
      if("wsum_gamma" %in% enrichment_method){
        print("wsum_gamma will start.")
        result_wsum_gamma=gamma_test(geneset=geneset,value=value)
        result_wsum_gamma$two=TRUE
        result_wsum_gamma$invert_A=ifelse(result_wsum_gamma$r<=0,TRUE,FALSE)
        result_wsum_gamma$invert_I=ifelse(result_wsum_gamma$r>0,TRUE,FALSE)
        result_wsum_gamma$P.Value.A=two2one(result_wsum_gamma$P.Value, two = as.logical(result_wsum_gamma$two), invert = as.logical(result_wsum_gamma$invert_A))
        result_wsum_gamma$P.Value.I=two2one(result_wsum_gamma$P.Value, two = as.logical(result_wsum_gamma$two), invert = as.logical(result_wsum_gamma$invert_I))
        result_wsum_gamma_A=result_wsum_gamma[,c("geneset_id","P.Value.A")]
        result_wsum_gamma_I=result_wsum_gamma[,c("geneset_id","P.Value.I")]
        colnames(result_wsum_gamma_A)[2]=paste0(vector[i],".wsum_gamma")
        colnames(result_wsum_gamma_I)[2]=paste0(vector[i],".wsum_gamma")
        rm(result_wsum_gamma)
      }
    }else {
      print("Not using wsum.")
    }
    ##########################################norm_wsum
    methods_to_check=c("norm_wsum_pearson","norm_wsum_kendall","norm_wsum_spearman","norm_wsum_lm","norm_wsum_biweight","norm_wsum_distance","norm_wsum_percentage","norm_wsum_blomqvist","norm_wsum_hoeffding","norm_wsum_gamma")
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
      if("norm_wsum_pearson" %in% enrichment_method){
        print("norm_wsum_pearson will start.")
        result_norm_wsum_pearson=pearson_test(geneset=geneset,value=value)
        result_norm_wsum_pearson$two=TRUE
        result_norm_wsum_pearson$invert_A=ifelse(result_norm_wsum_pearson$r<=0,TRUE,FALSE)
        result_norm_wsum_pearson$invert_I=ifelse(result_norm_wsum_pearson$r>0,TRUE,FALSE)
        result_norm_wsum_pearson$P.Value.A=two2one(result_norm_wsum_pearson$P.Value, two = as.logical(result_norm_wsum_pearson$two), invert = as.logical(result_norm_wsum_pearson$invert_A))
        result_norm_wsum_pearson$P.Value.I=two2one(result_norm_wsum_pearson$P.Value, two = as.logical(result_norm_wsum_pearson$two), invert = as.logical(result_norm_wsum_pearson$invert_I))
        result_norm_wsum_pearson_A=result_norm_wsum_pearson[,c("geneset_id","P.Value.A")]
        result_norm_wsum_pearson_I=result_norm_wsum_pearson[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wsum_pearson_A)[2]=paste0(vector[i],".norm_wsum_pearson")
        colnames(result_norm_wsum_pearson_I)[2]=paste0(vector[i],".norm_wsum_pearson")
        rm(result_norm_wsum_pearson)
      }
      if("norm_wsum_kendall" %in% enrichment_method){
        print("norm_wsum_kendall will start.")
        result_norm_wsum_kendall=kendall_test(geneset=geneset,value=value)
        result_norm_wsum_kendall$two=TRUE
        result_norm_wsum_kendall$invert_A=ifelse(result_norm_wsum_kendall$r<=0,TRUE,FALSE)
        result_norm_wsum_kendall$invert_I=ifelse(result_norm_wsum_kendall$r>0,TRUE,FALSE)
        result_norm_wsum_kendall$P.Value.A=two2one(result_norm_wsum_kendall$P.Value, two = as.logical(result_norm_wsum_kendall$two), invert = as.logical(result_norm_wsum_kendall$invert_A))
        result_norm_wsum_kendall$P.Value.I=two2one(result_norm_wsum_kendall$P.Value, two = as.logical(result_norm_wsum_kendall$two), invert = as.logical(result_norm_wsum_kendall$invert_I))
        result_norm_wsum_kendall_A=result_norm_wsum_kendall[,c("geneset_id","P.Value.A")]
        result_norm_wsum_kendall_I=result_norm_wsum_kendall[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wsum_kendall_A)[2]=paste0(vector[i],".norm_wsum_kendall")
        colnames(result_norm_wsum_kendall_I)[2]=paste0(vector[i],".norm_wsum_kendall")
        rm(result_norm_wsum_kendall)
      }
      if("norm_wsum_spearman" %in% enrichment_method){
        print("norm_wsum_spearman will start.")
        result_norm_wsum_spearman=spearman_test(geneset=geneset,value=value)
        result_norm_wsum_spearman$two=TRUE
        result_norm_wsum_spearman$invert_A=ifelse(result_norm_wsum_spearman$r<=0,TRUE,FALSE)
        result_norm_wsum_spearman$invert_I=ifelse(result_norm_wsum_spearman$r>0,TRUE,FALSE)
        result_norm_wsum_spearman$P.Value.A=two2one(result_norm_wsum_spearman$P.Value, two = as.logical(result_norm_wsum_spearman$two), invert = as.logical(result_norm_wsum_spearman$invert_A))
        result_norm_wsum_spearman$P.Value.I=two2one(result_norm_wsum_spearman$P.Value, two = as.logical(result_norm_wsum_spearman$two), invert = as.logical(result_norm_wsum_spearman$invert_I))
        result_norm_wsum_spearman_A=result_norm_wsum_spearman[,c("geneset_id","P.Value.A")]
        result_norm_wsum_spearman_I=result_norm_wsum_spearman[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wsum_spearman_A)[2]=paste0(vector[i],".norm_wsum_spearman")
        colnames(result_norm_wsum_spearman_I)[2]=paste0(vector[i],".norm_wsum_spearman")
        rm(result_norm_wsum_spearman)
      }
      if("norm_wsum_lm" %in% enrichment_method){
        print("norm_wsum_lm will start.")
        result_norm_wsum_lm=lm_test(geneset=geneset,value=value)
        result_norm_wsum_lm$two=TRUE
        result_norm_wsum_lm$invert_A=ifelse(result_norm_wsum_lm$r<=0,TRUE,FALSE)
        result_norm_wsum_lm$invert_I=ifelse(result_norm_wsum_lm$r>0,TRUE,FALSE)
        result_norm_wsum_lm$P.Value.A=two2one(result_norm_wsum_lm$P.Value, two = as.logical(result_norm_wsum_lm$two), invert = as.logical(result_norm_wsum_lm$invert_A))
        result_norm_wsum_lm$P.Value.I=two2one(result_norm_wsum_lm$P.Value, two = as.logical(result_norm_wsum_lm$two), invert = as.logical(result_norm_wsum_lm$invert_I))
        result_norm_wsum_lm_A=result_norm_wsum_lm[,c("geneset_id","P.Value.A")]
        result_norm_wsum_lm_I=result_norm_wsum_lm[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wsum_lm_A)[2]=paste0(vector[i],".norm_wsum_lm")
        colnames(result_norm_wsum_lm_I)[2]=paste0(vector[i],".norm_wsum_lm")
        rm(result_norm_wsum_lm)
      }
      if("norm_wsum_biweight" %in% enrichment_method){
        print("norm_wsum_biweight will start.")
        result_norm_wsum_biweight=biweight_test(geneset=geneset,value=value)
        result_norm_wsum_biweight$two=TRUE
        result_norm_wsum_biweight$invert_A=ifelse(result_norm_wsum_biweight$r<=0,TRUE,FALSE)
        result_norm_wsum_biweight$invert_I=ifelse(result_norm_wsum_biweight$r>0,TRUE,FALSE)
        result_norm_wsum_biweight$P.Value.A=two2one(result_norm_wsum_biweight$P.Value, two = as.logical(result_norm_wsum_biweight$two), invert = as.logical(result_norm_wsum_biweight$invert_A))
        result_norm_wsum_biweight$P.Value.I=two2one(result_norm_wsum_biweight$P.Value, two = as.logical(result_norm_wsum_biweight$two), invert = as.logical(result_norm_wsum_biweight$invert_I))
        result_norm_wsum_biweight_A=result_norm_wsum_biweight[,c("geneset_id","P.Value.A")]
        result_norm_wsum_biweight_I=result_norm_wsum_biweight[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wsum_biweight_A)[2]=paste0(vector[i],".norm_wsum_biweight")
        colnames(result_norm_wsum_biweight_I)[2]=paste0(vector[i],".norm_wsum_biweight")
        rm(result_norm_wsum_biweight)
      }
      if("norm_wsum_distance" %in% enrichment_method){
        print("norm_wsum_distance will start.")
        result_norm_wsum_distance=distance_test(geneset=geneset,value=value)
        result_norm_wsum_distance$two=TRUE
        result_norm_wsum_distance$invert_A=ifelse(result_norm_wsum_distance$r<=0,TRUE,FALSE)
        result_norm_wsum_distance$invert_I=ifelse(result_norm_wsum_distance$r>0,TRUE,FALSE)
        result_norm_wsum_distance$P.Value.A=two2one(result_norm_wsum_distance$P.Value, two = as.logical(result_norm_wsum_distance$two), invert = as.logical(result_norm_wsum_distance$invert_A))
        result_norm_wsum_distance$P.Value.I=two2one(result_norm_wsum_distance$P.Value, two = as.logical(result_norm_wsum_distance$two), invert = as.logical(result_norm_wsum_distance$invert_I))
        result_norm_wsum_distance_A=result_norm_wsum_distance[,c("geneset_id","P.Value.A")]
        result_norm_wsum_distance_I=result_norm_wsum_distance[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wsum_distance_A)[2]=paste0(vector[i],".norm_wsum_distance")
        colnames(result_norm_wsum_distance_I)[2]=paste0(vector[i],".norm_wsum_distance")
        rm(result_norm_wsum_distance)
      }
      if("norm_wsum_percentage" %in% enrichment_method){
        print("norm_wsum_percentage will start.")
        result_norm_wsum_percentage=percentage_test(geneset=geneset,value=value)
        result_norm_wsum_percentage$two=TRUE
        result_norm_wsum_percentage$invert_A=ifelse(result_norm_wsum_percentage$r<=0,TRUE,FALSE)
        result_norm_wsum_percentage$invert_I=ifelse(result_norm_wsum_percentage$r>0,TRUE,FALSE)
        result_norm_wsum_percentage$P.Value.A=two2one(result_norm_wsum_percentage$P.Value, two = as.logical(result_norm_wsum_percentage$two), invert = as.logical(result_norm_wsum_percentage$invert_A))
        result_norm_wsum_percentage$P.Value.I=two2one(result_norm_wsum_percentage$P.Value, two = as.logical(result_norm_wsum_percentage$two), invert = as.logical(result_norm_wsum_percentage$invert_I))
        result_norm_wsum_percentage_A=result_norm_wsum_percentage[,c("geneset_id","P.Value.A")]
        result_norm_wsum_percentage_I=result_norm_wsum_percentage[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wsum_percentage_A)[2]=paste0(vector[i],".norm_wsum_percentage")
        colnames(result_norm_wsum_percentage_I)[2]=paste0(vector[i],".norm_wsum_percentage")
        rm(result_norm_wsum_percentage)
      }
      if("norm_wsum_blomqvist" %in% enrichment_method){
        print("norm_wsum_blomqvist will start.")
        result_norm_wsum_blomqvist=blomqvist_test(geneset=geneset,value=value)
        result_norm_wsum_blomqvist$two=TRUE
        result_norm_wsum_blomqvist$invert_A=ifelse(result_norm_wsum_blomqvist$r<=0,TRUE,FALSE)
        result_norm_wsum_blomqvist$invert_I=ifelse(result_norm_wsum_blomqvist$r>0,TRUE,FALSE)
        result_norm_wsum_blomqvist$P.Value.A=two2one(result_norm_wsum_blomqvist$P.Value, two = as.logical(result_norm_wsum_blomqvist$two), invert = as.logical(result_norm_wsum_blomqvist$invert_A))
        result_norm_wsum_blomqvist$P.Value.I=two2one(result_norm_wsum_blomqvist$P.Value, two = as.logical(result_norm_wsum_blomqvist$two), invert = as.logical(result_norm_wsum_blomqvist$invert_I))
        result_norm_wsum_blomqvist_A=result_norm_wsum_blomqvist[,c("geneset_id","P.Value.A")]
        result_norm_wsum_blomqvist_I=result_norm_wsum_blomqvist[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wsum_blomqvist_A)[2]=paste0(vector[i],".norm_wsum_blomqvist")
        colnames(result_norm_wsum_blomqvist_I)[2]=paste0(vector[i],".norm_wsum_blomqvist")
        rm(result_norm_wsum_blomqvist)
      }
      if("norm_wsum_hoeffding" %in% enrichment_method){
        print("norm_wsum_hoeffding will start.")
        result_norm_wsum_hoeffding=hoeffding_test(geneset=geneset,value=value)
        result_norm_wsum_hoeffding$two=TRUE
        result_norm_wsum_hoeffding$invert_A=ifelse(result_norm_wsum_hoeffding$r<=0,TRUE,FALSE)
        result_norm_wsum_hoeffding$invert_I=ifelse(result_norm_wsum_hoeffding$r>0,TRUE,FALSE)
        result_norm_wsum_hoeffding$P.Value.A=two2one(result_norm_wsum_hoeffding$P.Value, two = as.logical(result_norm_wsum_hoeffding$two), invert = as.logical(result_norm_wsum_hoeffding$invert_A))
        result_norm_wsum_hoeffding$P.Value.I=two2one(result_norm_wsum_hoeffding$P.Value, two = as.logical(result_norm_wsum_hoeffding$two), invert = as.logical(result_norm_wsum_hoeffding$invert_I))
        result_norm_wsum_hoeffding_A=result_norm_wsum_hoeffding[,c("geneset_id","P.Value.A")]
        result_norm_wsum_hoeffding_I=result_norm_wsum_hoeffding[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wsum_hoeffding_A)[2]=paste0(vector[i],".norm_wsum_hoeffding")
        colnames(result_norm_wsum_hoeffding_I)[2]=paste0(vector[i],".norm_wsum_hoeffding")
        rm(result_norm_wsum_hoeffding)
      }
      if("norm_wsum_gamma" %in% enrichment_method){
        print("norm_wsum_gamma will start.")
        result_norm_wsum_gamma=gamma_test(geneset=geneset,value=value)
        result_norm_wsum_gamma$two=TRUE
        result_norm_wsum_gamma$invert_A=ifelse(result_norm_wsum_gamma$r<=0,TRUE,FALSE)
        result_norm_wsum_gamma$invert_I=ifelse(result_norm_wsum_gamma$r>0,TRUE,FALSE)
        result_norm_wsum_gamma$P.Value.A=two2one(result_norm_wsum_gamma$P.Value, two = as.logical(result_norm_wsum_gamma$two), invert = as.logical(result_norm_wsum_gamma$invert_A))
        result_norm_wsum_gamma$P.Value.I=two2one(result_norm_wsum_gamma$P.Value, two = as.logical(result_norm_wsum_gamma$two), invert = as.logical(result_norm_wsum_gamma$invert_I))
        result_norm_wsum_gamma_A=result_norm_wsum_gamma[,c("geneset_id","P.Value.A")]
        result_norm_wsum_gamma_I=result_norm_wsum_gamma[,c("geneset_id","P.Value.I")]
        colnames(result_norm_wsum_gamma_A)[2]=paste0(vector[i],".norm_wsum_gamma")
        colnames(result_norm_wsum_gamma_I)[2]=paste0(vector[i],".norm_wsum_gamma")
        rm(result_norm_wsum_gamma)
      }
    }else {
      print("Not using norm_wsum.")
    }
    ##########################################corr_wsum
    methods_to_check=c("corr_wsum_pearson","corr_wsum_kendall","corr_wsum_spearman","corr_wsum_lm","corr_wsum_biweight","corr_wsum_distance","corr_wsum_percentage","corr_wsum_blomqvist","corr_wsum_hoeffding","corr_wsum_gamma")
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
      if("corr_wsum_pearson" %in% enrichment_method){
        print("corr_wsum_pearson will start.")
        result_corr_wsum_pearson=pearson_test(geneset=geneset,value=value)
        result_corr_wsum_pearson$two=TRUE
        result_corr_wsum_pearson$invert_A=ifelse(result_corr_wsum_pearson$r<=0,TRUE,FALSE)
        result_corr_wsum_pearson$invert_I=ifelse(result_corr_wsum_pearson$r>0,TRUE,FALSE)
        result_corr_wsum_pearson$P.Value.A=two2one(result_corr_wsum_pearson$P.Value, two = as.logical(result_corr_wsum_pearson$two), invert = as.logical(result_corr_wsum_pearson$invert_A))
        result_corr_wsum_pearson$P.Value.I=two2one(result_corr_wsum_pearson$P.Value, two = as.logical(result_corr_wsum_pearson$two), invert = as.logical(result_corr_wsum_pearson$invert_I))
        result_corr_wsum_pearson_A=result_corr_wsum_pearson[,c("geneset_id","P.Value.A")]
        result_corr_wsum_pearson_I=result_corr_wsum_pearson[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wsum_pearson_A)[2]=paste0(vector[i],".corr_wsum_pearson")
        colnames(result_corr_wsum_pearson_I)[2]=paste0(vector[i],".corr_wsum_pearson")
        rm(result_corr_wsum_pearson)
      }
      if("corr_wsum_kendall" %in% enrichment_method){
        print("corr_wsum_kendall will start.")
        result_corr_wsum_kendall=kendall_test(geneset=geneset,value=value)
        result_corr_wsum_kendall$two=TRUE
        result_corr_wsum_kendall$invert_A=ifelse(result_corr_wsum_kendall$r<=0,TRUE,FALSE)
        result_corr_wsum_kendall$invert_I=ifelse(result_corr_wsum_kendall$r>0,TRUE,FALSE)
        result_corr_wsum_kendall$P.Value.A=two2one(result_corr_wsum_kendall$P.Value, two = as.logical(result_corr_wsum_kendall$two), invert = as.logical(result_corr_wsum_kendall$invert_A))
        result_corr_wsum_kendall$P.Value.I=two2one(result_corr_wsum_kendall$P.Value, two = as.logical(result_corr_wsum_kendall$two), invert = as.logical(result_corr_wsum_kendall$invert_I))
        result_corr_wsum_kendall_A=result_corr_wsum_kendall[,c("geneset_id","P.Value.A")]
        result_corr_wsum_kendall_I=result_corr_wsum_kendall[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wsum_kendall_A)[2]=paste0(vector[i],".corr_wsum_kendall")
        colnames(result_corr_wsum_kendall_I)[2]=paste0(vector[i],".corr_wsum_kendall")
        rm(result_corr_wsum_kendall)
      }
      if("corr_wsum_spearman" %in% enrichment_method){
        print("corr_wsum_spearman will start.")
        result_corr_wsum_spearman=spearman_test(geneset=geneset,value=value)
        result_corr_wsum_spearman$two=TRUE
        result_corr_wsum_spearman$invert_A=ifelse(result_corr_wsum_spearman$r<=0,TRUE,FALSE)
        result_corr_wsum_spearman$invert_I=ifelse(result_corr_wsum_spearman$r>0,TRUE,FALSE)
        result_corr_wsum_spearman$P.Value.A=two2one(result_corr_wsum_spearman$P.Value, two = as.logical(result_corr_wsum_spearman$two), invert = as.logical(result_corr_wsum_spearman$invert_A))
        result_corr_wsum_spearman$P.Value.I=two2one(result_corr_wsum_spearman$P.Value, two = as.logical(result_corr_wsum_spearman$two), invert = as.logical(result_corr_wsum_spearman$invert_I))
        result_corr_wsum_spearman_A=result_corr_wsum_spearman[,c("geneset_id","P.Value.A")]
        result_corr_wsum_spearman_I=result_corr_wsum_spearman[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wsum_spearman_A)[2]=paste0(vector[i],".corr_wsum_spearman")
        colnames(result_corr_wsum_spearman_I)[2]=paste0(vector[i],".corr_wsum_spearman")
        rm(result_corr_wsum_spearman)
      }
      if("corr_wsum_lm" %in% enrichment_method){
        print("corr_wsum_lm will start.")
        result_corr_wsum_lm=lm_test(geneset=geneset,value=value)
        result_corr_wsum_lm$two=TRUE
        result_corr_wsum_lm$invert_A=ifelse(result_corr_wsum_lm$r<=0,TRUE,FALSE)
        result_corr_wsum_lm$invert_I=ifelse(result_corr_wsum_lm$r>0,TRUE,FALSE)
        result_corr_wsum_lm$P.Value.A=two2one(result_corr_wsum_lm$P.Value, two = as.logical(result_corr_wsum_lm$two), invert = as.logical(result_corr_wsum_lm$invert_A))
        result_corr_wsum_lm$P.Value.I=two2one(result_corr_wsum_lm$P.Value, two = as.logical(result_corr_wsum_lm$two), invert = as.logical(result_corr_wsum_lm$invert_I))
        result_corr_wsum_lm_A=result_corr_wsum_lm[,c("geneset_id","P.Value.A")]
        result_corr_wsum_lm_I=result_corr_wsum_lm[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wsum_lm_A)[2]=paste0(vector[i],".corr_wsum_lm")
        colnames(result_corr_wsum_lm_I)[2]=paste0(vector[i],".corr_wsum_lm")
        rm(result_corr_wsum_lm)
      }
      if("corr_wsum_biweight" %in% enrichment_method){
        print("corr_wsum_biweight will start.")
        result_corr_wsum_biweight=biweight_test(geneset=geneset,value=value)
        result_corr_wsum_biweight$two=TRUE
        result_corr_wsum_biweight$invert_A=ifelse(result_corr_wsum_biweight$r<=0,TRUE,FALSE)
        result_corr_wsum_biweight$invert_I=ifelse(result_corr_wsum_biweight$r>0,TRUE,FALSE)
        result_corr_wsum_biweight$P.Value.A=two2one(result_corr_wsum_biweight$P.Value, two = as.logical(result_corr_wsum_biweight$two), invert = as.logical(result_corr_wsum_biweight$invert_A))
        result_corr_wsum_biweight$P.Value.I=two2one(result_corr_wsum_biweight$P.Value, two = as.logical(result_corr_wsum_biweight$two), invert = as.logical(result_corr_wsum_biweight$invert_I))
        result_corr_wsum_biweight_A=result_corr_wsum_biweight[,c("geneset_id","P.Value.A")]
        result_corr_wsum_biweight_I=result_corr_wsum_biweight[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wsum_biweight_A)[2]=paste0(vector[i],".corr_wsum_biweight")
        colnames(result_corr_wsum_biweight_I)[2]=paste0(vector[i],".corr_wsum_biweight")
        rm(result_corr_wsum_biweight)
      }
      if("corr_wsum_distance" %in% enrichment_method){
        print("corr_wsum_distance will start.")
        result_corr_wsum_distance=distance_test(geneset=geneset,value=value)
        result_corr_wsum_distance$two=TRUE
        result_corr_wsum_distance$invert_A=ifelse(result_corr_wsum_distance$r<=0,TRUE,FALSE)
        result_corr_wsum_distance$invert_I=ifelse(result_corr_wsum_distance$r>0,TRUE,FALSE)
        result_corr_wsum_distance$P.Value.A=two2one(result_corr_wsum_distance$P.Value, two = as.logical(result_corr_wsum_distance$two), invert = as.logical(result_corr_wsum_distance$invert_A))
        result_corr_wsum_distance$P.Value.I=two2one(result_corr_wsum_distance$P.Value, two = as.logical(result_corr_wsum_distance$two), invert = as.logical(result_corr_wsum_distance$invert_I))
        result_corr_wsum_distance_A=result_corr_wsum_distance[,c("geneset_id","P.Value.A")]
        result_corr_wsum_distance_I=result_corr_wsum_distance[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wsum_distance_A)[2]=paste0(vector[i],".corr_wsum_distance")
        colnames(result_corr_wsum_distance_I)[2]=paste0(vector[i],".corr_wsum_distance")
        rm(result_corr_wsum_distance)
      }
      if("corr_wsum_percentage" %in% enrichment_method){
        print("corr_wsum_percentage will start.")
        result_corr_wsum_percentage=percentage_test(geneset=geneset,value=value)
        result_corr_wsum_percentage$two=TRUE
        result_corr_wsum_percentage$invert_A=ifelse(result_corr_wsum_percentage$r<=0,TRUE,FALSE)
        result_corr_wsum_percentage$invert_I=ifelse(result_corr_wsum_percentage$r>0,TRUE,FALSE)
        result_corr_wsum_percentage$P.Value.A=two2one(result_corr_wsum_percentage$P.Value, two = as.logical(result_corr_wsum_percentage$two), invert = as.logical(result_corr_wsum_percentage$invert_A))
        result_corr_wsum_percentage$P.Value.I=two2one(result_corr_wsum_percentage$P.Value, two = as.logical(result_corr_wsum_percentage$two), invert = as.logical(result_corr_wsum_percentage$invert_I))
        result_corr_wsum_percentage_A=result_corr_wsum_percentage[,c("geneset_id","P.Value.A")]
        result_corr_wsum_percentage_I=result_corr_wsum_percentage[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wsum_percentage_A)[2]=paste0(vector[i],".corr_wsum_percentage")
        colnames(result_corr_wsum_percentage_I)[2]=paste0(vector[i],".corr_wsum_percentage")
        rm(result_corr_wsum_percentage)
      }
      if("corr_wsum_blomqvist" %in% enrichment_method){
        print("corr_wsum_blomqvist will start.")
        result_corr_wsum_blomqvist=blomqvist_test(geneset=geneset,value=value)
        result_corr_wsum_blomqvist$two=TRUE
        result_corr_wsum_blomqvist$invert_A=ifelse(result_corr_wsum_blomqvist$r<=0,TRUE,FALSE)
        result_corr_wsum_blomqvist$invert_I=ifelse(result_corr_wsum_blomqvist$r>0,TRUE,FALSE)
        result_corr_wsum_blomqvist$P.Value.A=two2one(result_corr_wsum_blomqvist$P.Value, two = as.logical(result_corr_wsum_blomqvist$two), invert = as.logical(result_corr_wsum_blomqvist$invert_A))
        result_corr_wsum_blomqvist$P.Value.I=two2one(result_corr_wsum_blomqvist$P.Value, two = as.logical(result_corr_wsum_blomqvist$two), invert = as.logical(result_corr_wsum_blomqvist$invert_I))
        result_corr_wsum_blomqvist_A=result_corr_wsum_blomqvist[,c("geneset_id","P.Value.A")]
        result_corr_wsum_blomqvist_I=result_corr_wsum_blomqvist[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wsum_blomqvist_A)[2]=paste0(vector[i],".corr_wsum_blomqvist")
        colnames(result_corr_wsum_blomqvist_I)[2]=paste0(vector[i],".corr_wsum_blomqvist")
        rm(result_corr_wsum_blomqvist)
      }
      if("corr_wsum_hoeffding" %in% enrichment_method){
        print("corr_wsum_hoeffding will start.")
        result_corr_wsum_hoeffding=hoeffding_test(geneset=geneset,value=value)
        result_corr_wsum_hoeffding$two=TRUE
        result_corr_wsum_hoeffding$invert_A=ifelse(result_corr_wsum_hoeffding$r<=0,TRUE,FALSE)
        result_corr_wsum_hoeffding$invert_I=ifelse(result_corr_wsum_hoeffding$r>0,TRUE,FALSE)
        result_corr_wsum_hoeffding$P.Value.A=two2one(result_corr_wsum_hoeffding$P.Value, two = as.logical(result_corr_wsum_hoeffding$two), invert = as.logical(result_corr_wsum_hoeffding$invert_A))
        result_corr_wsum_hoeffding$P.Value.I=two2one(result_corr_wsum_hoeffding$P.Value, two = as.logical(result_corr_wsum_hoeffding$two), invert = as.logical(result_corr_wsum_hoeffding$invert_I))
        result_corr_wsum_hoeffding_A=result_corr_wsum_hoeffding[,c("geneset_id","P.Value.A")]
        result_corr_wsum_hoeffding_I=result_corr_wsum_hoeffding[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wsum_hoeffding_A)[2]=paste0(vector[i],".corr_wsum_hoeffding")
        colnames(result_corr_wsum_hoeffding_I)[2]=paste0(vector[i],".corr_wsum_hoeffding")
        rm(result_corr_wsum_hoeffding)
      }
      if("corr_wsum_gamma" %in% enrichment_method){
        print("corr_wsum_gamma will start.")
        result_corr_wsum_gamma=gamma_test(geneset=geneset,value=value)
        result_corr_wsum_gamma$two=TRUE
        result_corr_wsum_gamma$invert_A=ifelse(result_corr_wsum_gamma$r<=0,TRUE,FALSE)
        result_corr_wsum_gamma$invert_I=ifelse(result_corr_wsum_gamma$r>0,TRUE,FALSE)
        result_corr_wsum_gamma$P.Value.A=two2one(result_corr_wsum_gamma$P.Value, two = as.logical(result_corr_wsum_gamma$two), invert = as.logical(result_corr_wsum_gamma$invert_A))
        result_corr_wsum_gamma$P.Value.I=two2one(result_corr_wsum_gamma$P.Value, two = as.logical(result_corr_wsum_gamma$two), invert = as.logical(result_corr_wsum_gamma$invert_I))
        result_corr_wsum_gamma_A=result_corr_wsum_gamma[,c("geneset_id","P.Value.A")]
        result_corr_wsum_gamma_I=result_corr_wsum_gamma[,c("geneset_id","P.Value.I")]
        colnames(result_corr_wsum_gamma_A)[2]=paste0(vector[i],".corr_wsum_gamma")
        colnames(result_corr_wsum_gamma_I)[2]=paste0(vector[i],".corr_wsum_gamma")
        rm(result_corr_wsum_gamma)
      }
    }else {
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
