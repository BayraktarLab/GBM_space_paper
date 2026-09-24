#-------------------------------------------------
# Shearwater-like filter for WGS (post-CaVEMan)
# Tim Coorens - November 2018
#-------------------------------------------------
options(stringsAsFactors=F)

#-------------------------------------------------
# Libraries
#-------------------------------------------------

library("GenomicRanges")
library("deepSNV")
library("Rsamtools")

logbb = deepSNV:::logbb
dbetabinom = VGAM::dbetabinom

#-------------------------------------------------
# Functions
#-------------------------------------------------

estimateRho_gridml = function(x, mu) {
  # Estimate rho by MLE grid approach
  rhovec = 10^seq(-6,-0.5,by=0.05) # rho will be bounded within 1e-6 and 0.32
  mm = x[,2]
  #cov = c(x[,3:4])+c(x[,1:2])
  cov = c(x[,1])
  ll = sapply(rhovec, function(rhoj) sum(dbetabinom(x=mm, size=cov, rho=rhoj, prob=mu, log=T)))
  rhovec[ll==max(ll)][1]
}
#path_prefix="/lustre/scratch126/casm/team274sb/aw35/bilateral_WT"
shearwater_probability=function(patient, save=NULL,allVar,all_counts,norm_all_counts,case_samples, rho=10^-3){
  # Function to calculcate probability of presence of mutation based on Shearwarer
  # 'patient' is the name of the patient-specific subfolder
  # 'path_prefix' is any prefix to the path necessary to find that subfolder
  # 'rho' is the constant for the overdispersion parameter. If rho=NULL, calculate
  # it from the data (much slower)
  # 'save' is path for output. If NULL, returns matrix
  # output of this function will be a matrix (muts by samples) of p values

  
  # A file of mutations in patient subdirectory (format: Chr_Ref_Pos_Alt)
  #Muts_patient = read.table(paste0(path_prefix,"/mutations/germline_flt/",patient,"_mut_germline_flt.txt"), header=T)#[,1]
  #Muts_patient=paste(allVar$Chr,allVar$Pos, allVar$Ref, allVar$Alt, sep="_")
  # A file with the sample names belonging to this patient
  #all_samples=read.table(paste0(path_prefix,"/samples_bilateral_WT.txt"), header=T, sep='\t')[,1]
  #case_samples=all_samples[grepl(patient, all_samples)]
  
  # Set up pval matrix
  pval_mat = matrix(1,ncol=nrow(all_counts),nrow=ncol(all_counts))
  rownames(pval_mat)=colnames(all_counts)
  colnames(pval_mat)=rownames(all_counts)
  
  coords_proj = colnames(all_counts)
  Alt=allVar$Alt
  Ref=allVar$Ref
  
  for (s in case_samples){
    rho_est=rep(NA,length(coords_proj))
    test_counts = all_counts[s,coords_proj,]
    for (k in 1:length(coords_proj)) {
      n = sum(test_counts[coords_proj[k],])
      x = test_counts[coords_proj[k],Alt[k]]
      
      N_indiv = rowSums(norm_all_counts[,coords_proj[k],])
      X_indiv = norm_all_counts[,coords_proj[k],c("A","C","G","T")!=Ref[k]]
      pseudo = .Machine$double.eps    
      N=sum(N_indiv)
      X=sum(X_indiv)
      
      mu = max(X,pseudo)/max(N,pseudo)
      counts = cbind(N,X)
      if(is.null(rho)) rho = estimateRho_gridml(counts,mu)
      rdisp = (1 - rho)/rho
      
      prob0 = (X + x)/(N + n); prob0[prob0==0] = pseudo
      prob1s = x/(n+pseudo); prob1s[prob1s==0] = pseudo
      prob1c = X/(N+pseudo); prob1c[prob1c==0] = pseudo
      
      prob1s = pmax(prob1s,prob1c) # Min error rate is that of the population (one-sided test)
      nu0 = prob0 * rdisp; nu1s = prob1s * rdisp; nu1c = prob1c * rdisp; 
      
      # Likelihood-Ratio Tests
      LL = logbb(x, n, nu0, rdisp) + logbb(X, N, nu0, rdisp) - logbb(x, n, nu1s, rdisp) - logbb(X, N, nu1c, rdisp)
      pvals = pchisq(-2*LL, df=1, lower.tail=F)/2 # We divide by 2 as we are performing a 1-sided test
      # Saving the result
      pval_mat[k,s] = pvals
    } 
  }
  if(is.null(save)){
    return(pval_mat)
  }else{
    write.table(pval_mat,save)
    return(pval_mat)
  }
}





#-------------------------------------------------
# Shearwater-like filter for WGS (post-CaVEMan)
# Wrapper functions - Mi 2023
#-------------------------------------------------

extract_alleleCount = function(outDir,cgpVAF_output_fp,
                               alleleCount_outDir,samples_PDID){
  
  ## check for proper inputs
  if(is.null(samples_PDID)){
    stop('Please provide samples_PDID')
  }
  if(is.null(cgpVAF_output_fp) | is.null(alleleCount_outDir)){
    stop('Please provide cgpVAF_output_fp and alleleCount_outDir')
  }
  if(!file.exists(cgpVAF_output_fp)){
    stop('cgpVAF_output_fp does not exist. Please check!')
  }
  if(!dir.exists(alleleCount_outDir)){
    message('Creating alleleCount_outDir')
    dir.create(alleleCount_outDir,recursive = T)
  }
  
  
  ##-------------------------------------------------------------------------------##
  ##    STEP2: alleleCount at each variant position across panel of normal samples ##
  ##-------------------------------------------------------------------------------##
  
  print('STEP2: alleleCount at each variant position across panel of normal samples')
  
  # Read in cgpVAF output 
  cgpVAF_output = read.table(cgpVAF_output_fp,header = TRUE, sep = '\t', stringsAsFactors = F, comment.char = "",skip = 64)
  colnames(cgpVAF_output)[which(colnames(cgpVAF_output)=="VariantID")]="ID"
  row.names(cgpVAF_output)=paste0(sub('.*chr', '', cgpVAF_output$Chrom), "_", cgpVAF_output$Pos)
  
  
  
  allele_count = cgpVAF_output[,c('Chrom','Pos')]
  
  # Extract base count for each sample
  for (i in seq_along(samples_PDID)) {
    sample = samples_PDID[i]
    print(sample)
    # Check if it is present in cgpVAF output
    if(sum(grepl(paste0(sample,"_"),colnames(cgpVAF_output))) == 0){
      warning(sprintf('No cgpVAF output for sample %s',sample))
      next()
    }else{
      mut_temp=cgpVAF_output[,grep(paste0(sample,"_"), colnames(cgpVAF_output))]  
      #reset the allele counts to NA before starting a new sample. Not strictly necessary but doing it for my sanity/piece of mind
      allele_count$Count_A = mut_temp[,grepl('FAZ',colnames(mut_temp))] + mut_temp[,grepl('RAZ',colnames(mut_temp))]
      allele_count$Count_C = mut_temp[,grepl('FCZ',colnames(mut_temp))] + mut_temp[,grepl('RCZ',colnames(mut_temp))]
      allele_count$Count_G = mut_temp[,grepl('FGZ',colnames(mut_temp))] + mut_temp[,grepl('RGZ',colnames(mut_temp))]
      allele_count$Count_T = mut_temp[,grepl('FTZ',colnames(mut_temp))] + mut_temp[,grepl('RTZ',colnames(mut_temp))]
    }
    
    write.table(allele_count, file = file.path(alleleCount_outDir,paste0(sample, "_alleleCounts.txt")),
                sep = "\t", row.names = F, quote = F, col.names = T)
  }
  
  print('STEP2: Completed')
}









runShearWaterLike = function(case_samples_PDID,normal_samples_PDID,
                             outDir,variants,alleleCount_outDir,
                             FDR=0.001,patient=NULL){
  
  ## check for proper inputs
  if(is.null(case_samples_PDID)){
    stop('Please provide case_samples_PDID')
  }
  if(is.null(normal_samples_PDID)){
    stop('Please provide normal_samples_PDID')
  }
  
  
  ##----------------------------------------------------------------------##
  ##    STEP3: Run Tim's shearwater-like filter to remove false mutations ##
  ##----------------------------------------------------------------------##
  
  print('STEP3: Run Tims shearwater-like filter to remove false mutations')
  
  ##---- Inputs
  # Vector of all sample names
  samples_ID=case_samples_PDID
  
  # Vector of normal panel (only blood from other project [3110])
  normal_samples_PDID = normal_samples_PDID
  
  # Table of all variants to be considered
  colnames(variants) = c('Chr','Pos','Ref','Alt')
  coords=paste(variants$Chr,variants$Pos,sep="_")
  rownames(variants) = coords
  
  # Read in data from AlleleCounter/cgpVAF
  all_counts = array(0,dim=c(length(samples_ID),length(coords),4),
                     dimnames=list(samples_ID,coords,c("A","C","G","T")))
  print(length(samples_ID))
  for (k in 1:length(samples_ID)){
    #Read in allele counts per sample
    if(file.exists(file.path(alleleCount_outDir, paste0(samples_ID[k],"_alleleCounts.txt")))){
      print(samples_ID[k])
      data=read.table(file.path(alleleCount_outDir, paste0(samples_ID[k],"_alleleCounts.txt")),comment.char = '',header=T)
      rownames(data)=paste(data$Chrom,data$Pos,sep="_")
      all_counts[k,,]=as.matrix(data[coords,3:6])
    }
  }
  
  
  # Read in data from AlleleCounter/cgpVAF for all samples in the normal panel
  norm_all_counts = array(0,dim=c(length(normal_samples_PDID),length(coords),4),
                          dimnames=list(normal_samples_PDID,coords,c("A","C","G","T")))
  
  for (k in seq_along(normal_samples_PDID)){
    sample = normal_samples_PDID[k]
    #Read in allele counts per sample
    if(file.exists(file.path(alleleCount_outDir, paste0(sample,"_alleleCounts.txt")))){
      print(normal_samples_PDID[k])
      data=read.table(file.path(alleleCount_outDir, paste0(sample,"_alleleCounts.txt")),comment.char = '',header=T)
      rownames(data)=paste(data$Chrom,data$Pos,sep="_")
      norm_all_counts[k,,]=as.matrix(data[coords,3:6])
    }
  }
  
  
  
  
  #------- Run shearwater filter
  if(!is.null(patient)){
    outFile = file.path(outDir,paste0(patient,"_shearwater_pval_mat.txt"))
  }else{
    outFile = file.path(outDir,"shearwater_pval_mat.txt")
  }
  pval_mat = shearwater_probability(patient, save=outFile,allVar=variants,
                                    all_counts=all_counts,norm_all_counts=norm_all_counts,case_samples=samples_ID, rho=10^-3)
  
  
  qval_mat=apply(pval_mat,2,function(x) p.adjust(x,method="BH",n = length(pval_mat)))
  # Note that if MTR=DEP, the p-value will be so low that it is listed as NA(n). Change these values to 0
  qval_mat[which(qval_mat=="NaN")]=0
  row.names(qval_mat)=sub('.*chr', '', row.names(qval_mat))
  
  # Add the qval info to the mutations object (make sure rows are in the same order first)
  mutations = variants
  row.names(mutations)=paste0(sub('.*chr', '', mutations$Chr), "_", mutations$Pos)#, "_", mutations$Ref,"_", mutations$Alt)
  table(row.names(mutations) %in% row.names(qval_mat))
  reorder_idx=match(row.names(mutations), row.names(qval_mat))
  qval_mat=qval_mat[reorder_idx, ]
  mutations=cbind(mutations, qval_mat)
  
  print('STEP3: Shearwater-like filter completed!')
  
  
  return(mutations)
}



