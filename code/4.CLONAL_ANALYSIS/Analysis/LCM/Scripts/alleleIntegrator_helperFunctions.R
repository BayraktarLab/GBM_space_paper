##--------------------------------------------##
##   Helper Functions for alleleIntegrator  ####
##--------------------------------------------##



#' @param btb.fp file path to btb output ".../xx.battenberg.subclones.txt.gz'
#' @param PDID PDID of the DNA sample being processed
#' @param tgtChrs chromosomes to consider
#' @param longFormat output data.frame with each segment represented by 2 rows, with Start and Stop position
#' @param minSegLen Minimum clonal CN segment length to consider. Any segments with shorter length will be ignored
#' @param subCl.minSegLen Minimum sub-clonal CN segment length to consider. Any segments with shorter length will be ignored
#' @param keepSubClonalSegs Return sub-clonal CN segments if TRUE
#' @param removeBalancedSegs Remove CN segments (clonal and sub-clonal) if major:minor allele ratio is 1:1 (i.e no allelic imbalances)
#' @param method Method to combine CN segments. "allelicRatio" = combine adjacent segments with the same allelic imbalances

processBTBsummary <- function(
    btb.fp,
    PDID,
    tgtChrs = c(1:23, "X"),
    longFormat = TRUE,
    minSegLen = 1e6,
    subCl.minSegLen = 1e7,
    keepSubClonalSegs = FALSE,
    removeBalancedSegs = FALSE,
    method = c("totalCN", "allelicRatio", "config"),
    chromInfo.fp = "/lustre/scratch125/cellgen/behjati/mt22/generalResources/chrom_abspos_kb.txt",
    summary.sep = ",",
    summary.header = FALSE,
    summary.colnames = c(
        "idx", "chr", "startpos", "endpos",
        "nMaj1_A", "nMin1_A", "nMaj2_A", "nMin2_A"
    )
) {
    method <- match.arg(method)

    if (!requireNamespace("GenomicRanges", quietly = TRUE) ||
        !requireNamespace("IRanges", quietly = TRUE)) {
        stop("Packages GenomicRanges and IRanges are required.")
    }
    if (longFormat && !requireNamespace("tidyr", quietly = TRUE)) {
        stop("Package tidyr is required when longFormat = TRUE.")
    }

    chromInfo <- utils::read.delim(chromInfo.fp, sep = "\t")

    # Read btb.fp itself.  Unlike the original function, do not replace
    # 'summary.csv' with 'subclones.txt.gz'.
    btb <- utils::read.table(
        btb.fp,
        sep = summary.sep,
        header = summary.header,
        stringsAsFactors = FALSE,
        na.strings = c("NA", "NaN", ".", ""),
        check.names = FALSE,
        comment.char = "",
        quote = "\""
    )

    if (!summary.header) {
        if (ncol(btb) != length(summary.colnames)) {
            stop(
                "Expected ", length(summary.colnames),
                " summary columns but found ", ncol(btb), "."
            )
        }
        names(btb) <- summary.colnames
    }

    required <- c(
        "chr", "startpos", "endpos",
        "nMaj1_A", "nMin1_A", "nMaj2_A", "nMin2_A"
    )
    missing.cols <- setdiff(required, names(btb))
    if (length(missing.cols) > 0L) {
        stop("Missing summary columns: ", paste(missing.cols, collapse = ", "))
    }
    if (!"idx" %in% names(btb)) {
        btb$idx <- seq_len(nrow(btb))
    }

    numeric.cols <- c(
        "idx", "startpos", "endpos",
        "nMaj1_A", "nMin1_A", "nMaj2_A", "nMin2_A"
    )
    btb[numeric.cols] <- lapply(btb[numeric.cols], function(x) {
        suppressWarnings(as.numeric(x))
    })

    normalise_chr <- function(x) {
        x <- sub("^chr", "", as.character(x), ignore.case = TRUE)
        x[toupper(x) == "X"] <- "23"
        suppressWarnings(as.integer(x))
    }

    btb$chr <- normalise_chr(btb$chr)
    target.chrs <- unique(normalise_chr(tgtChrs))
    target.chrs <- target.chrs[!is.na(target.chrs)]

    valid <- !is.na(btb$chr) &
        !is.na(btb$startpos) & !is.na(btb$endpos) &
        btb$startpos >= 1 & btb$endpos >= btb$startpos &
        !is.na(btb$nMaj1_A) & !is.na(btb$nMin1_A)
    if (any(!valid)) {
        warning(sum(!valid), " invalid summary rows were removed.")
        btb <- btb[valid, , drop = FALSE]
    }
    btb <- btb[btb$chr %in% target.chrs, , drop = FALSE]
    if (nrow(btb) == 0L) {
        stop("No valid summary segments were found on the target chromosomes.")
    }

    btb$posID <- paste0(btb$chr, ":", btb$startpos, "_", btb$endpos)
    has.state2 <- !is.na(btb$nMaj2_A) & !is.na(btb$nMin2_A)
    incomplete.state2 <- xor(is.na(btb$nMaj2_A), is.na(btb$nMin2_A))
    if (any(incomplete.state2)) {
        warning(
            sum(incomplete.state2),
            " rows have only one state-2 copy-number value; state 2 was ignored."
        )
    }

    safe_ratio <- function(major, minor) {
        total <- major + minor
        ans <- major / total
        ans[total == 0] <- 0
        ans
    }

    # Every summary row contributes its state-1 segment.  This is the key
    # difference from the old subclones.txt.gz logic: no cellular fractions
    # are invented, and state 1 is not discarded merely because state 2 exists.
    main.segs <- GenomicRanges::GRanges(
        seqnames = btb$chr,
        ranges = IRanges::IRanges(btb$startpos, btb$endpos),
        chr = btb$chr,
        sourceIdx = btb$idx,
        patNum = btb$nMin1_A,
        matNum = btb$nMaj1_A,
        totCN = btb$nMaj1_A + btb$nMin1_A,
        tumFrac = safe_ratio(btb$nMaj1_A, btb$nMin1_A),
        clonalType = ifelse(has.state2, "maj", "none")
    )

    sub.segs <- GenomicRanges::GRanges()
    if (keepSubClonalSegs && any(has.state2)) {
        s <- btb[has.state2, , drop = FALSE]
        sub.segs <- GenomicRanges::GRanges(
            seqnames = s$chr,
            ranges = IRanges::IRanges(s$startpos, s$endpos),
            chr = s$chr,
            sourceIdx = s$idx,
            patNum = s$nMin2_A,
            matNum = s$nMaj2_A,
            totCN = s$nMaj2_A + s$nMin2_A,
            tumFrac = safe_ratio(s$nMaj2_A, s$nMin2_A),
            clonalType = "sub"
        )
        sub.segs <- sub.segs[IRanges::width(sub.segs) >= subCl.minSegLen]
    }

    # Accumulate merged segments across chromosomes.  Separate accumulator
    # names avoid the original function's clonalSegs overwrite bug.
    merged.main <- GenomicRanges::GRanges()

    for (chr in target.chrs) {
        chr.main <- main.segs[GenomicRanges::seqnames(main.segs) == chr]
        if (length(chr.main) == 0L) {
            next
        }

        chr.segs <- chr.main

        # Retain the original convention: uncovered regions are diploid.
        chr.info <- chromInfo[chromInfo$chrom == chr, , drop = FALSE]
        if (nrow(chr.info) > 0L) {
            chrom.len <- max(chr.info$end, na.rm = TRUE) * 1000
            whole.chr <- GenomicRanges::GRanges(
                seqnames = chr,
                ranges = IRanges::IRanges(1, chrom.len)
            )
            covered <- GenomicRanges::reduce(chr.main, ignore.strand = TRUE)
            gap.segs <- GenomicRanges::setdiff(
                whole.chr, covered, ignore.strand = TRUE
            )
            if (length(gap.segs) > 0L) {
                S4Vectors::mcols(gap.segs) <- S4Vectors::DataFrame(
                    chr = rep(chr, length(gap.segs)),
                    sourceIdx = rep(NA_integer_, length(gap.segs)),
                    patNum = rep(1, length(gap.segs)),
                    matNum = rep(1, length(gap.segs)),
                    totCN = rep(2, length(gap.segs)),
                    tumFrac = rep(0.5, length(gap.segs)),
                    clonalType = rep("none", length(gap.segs))
                )
                chr.segs <- c(chr.segs, gap.segs)
            }
        }

        chr.segs$configuration <- paste0(
            chr.segs$totCN, ":", chr.segs$patNum
        )

        if (method == "totalCN") {
            group.value <- as.character(chr.segs$totCN)
        } else if (method == "allelicRatio") {
            group.value <- as.character(chr.segs$tumFrac)
        } else {
            group.value <- chr.segs$configuration
        }

        for (value in unique(group.value)) {
            z <- chr.segs[group.value == value]
            z.reduced <- GenomicRanges::reduce(z, ignore.strand = TRUE)

            if (method == "totalCN") {
                total <- unique(z$totCN)
                S4Vectors::mcols(z.reduced) <- S4Vectors::DataFrame(
                    chr = rep(chr, length(z.reduced)),
                    sourceIdx = rep(NA_integer_, length(z.reduced)),
                    patNum = rep(NA_real_, length(z.reduced)),
                    matNum = rep(NA_real_, length(z.reduced)),
                    totCN = rep(total[[1L]], length(z.reduced)),
                    tumFrac = rep(NA_real_, length(z.reduced)),
                    clonalType = rep("none", length(z.reduced))
                )
            } else if (method == "allelicRatio") {
                ratio <- unique(z$tumFrac)
                S4Vectors::mcols(z.reduced) <- S4Vectors::DataFrame(
                    chr = rep(chr, length(z.reduced)),
                    sourceIdx = rep(NA_integer_, length(z.reduced)),
                    patNum = rep(NA_real_, length(z.reduced)),
                    matNum = rep(NA_real_, length(z.reduced)),
                    totCN = rep(NA_real_, length(z.reduced)),
                    tumFrac = rep(ratio[[1L]], length(z.reduced)),
                    clonalType = rep("none", length(z.reduced))
                )
            } else {
                pieces <- strsplit(value, ":", fixed = TRUE)[[1L]]
                total <- as.numeric(pieces[[1L]])
                minor <- as.numeric(pieces[[2L]])
                S4Vectors::mcols(z.reduced) <- S4Vectors::DataFrame(
                    chr = rep(chr, length(z.reduced)),
                    sourceIdx = rep(NA_integer_, length(z.reduced)),
                    patNum = rep(minor, length(z.reduced)),
                    matNum = rep(total - minor, length(z.reduced)),
                    totCN = rep(total, length(z.reduced)),
                    tumFrac = rep(safe_ratio(total - minor, minor), length(z.reduced)),
                    clonalType = rep("none", length(z.reduced))
                )
            }

            merged.main <- c(merged.main, z.reduced)
        }
    }

    final.segs <- merged.main
    if (keepSubClonalSegs && length(sub.segs) > 0L) {
        final.segs <- c(final.segs, sub.segs)
    }

    final.segs <- final.segs[IRanges::width(final.segs) >= minSegLen]
    final.segs$Start <- IRanges::start(final.segs)
    final.segs$Stop <- IRanges::end(final.segs)
    final.segs$posID <- paste0(
        GenomicRanges::seqnames(final.segs), ":",
        final.segs$Start, "_", final.segs$Stop
    )

    if (removeBalancedSegs && length(final.segs) > 0L) {
        balanced <- !is.na(final.segs$tumFrac) & final.segs$tumFrac == 0.5
        final.segs <- final.segs[!balanced]
    }

    final.segs <- final.segs[order(final.segs$chr, final.segs$Start)]

    if (!longFormat) {
        final.segs$idx <- seq_along(final.segs)
        names(final.segs) <- as.character(GenomicRanges::seqnames(final.segs))
        return(final.segs)
    }

    data <- as.data.frame(S4Vectors::mcols(final.segs))
    out <- tidyr::pivot_longer(
        data,
        cols = c("Start", "Stop"),
        names_to = "posType",
        values_to = "pos"
    )
    out$abspos_kb <- out$pos / 1000

    for (r in seq_len(nrow(out))) {
        chr <- out$chr[[r]]
        if (chr > 1L) {
            prior <- chromInfo[
                chromInfo$chrom == chr - 1L & chromInfo$arm == "q",
                "abspos_kb"
            ]
            if (length(prior) > 0L) {
                out$abspos_kb[[r]] <- out$abspos_kb[[r]] + prior[[1L]]
            }
        }
    }

    out <- out[order(out$chr, out$abspos_kb), , drop = FALSE]
    out$idx <- seq_len(nrow(out))
    rownames(out) <- NULL
    out
}


processBTB = function(btb.fp,PDID,tgtChrs=c(1:23,'X'),longFormat=T,
                      minSegLen=1e6,subCl.minSegLen=1e7,keepSubClonalSegs=F,
                      removeBalancedSegs = FALSE,
                      method = c('totalCN','allelicRatio','config')){
  require(GenomicRanges)
  chromInfo = read.delim('/lustre/scratch125/cellgen/behjati/mt22/generalResources/chrom_abspos_kb.txt',sep = '\t')
  
  ## Read in CN profile for both major and minor clone from Battenberg ####
  btb = gsub('summary.csv','subclones.txt.gz',btb.fp)
  btb = read.delim(btb,sep = '\t')
  btb$posID = paste0(btb$chr,':',btb$startpos,'_',btb$endpos)
  btb$idx = c(1:nrow(btb))
  
  ## Remove clonality if subclone fraction is only <= 5%
  w = which(!is.na(btb$frac2_A) & btb$frac2_A < 0.05)
  btb[w,c("nMaj2_A","nMin2_A","frac2_A")] = NA
  if(sum(btb$frac1_A < 0.05) > 0){
    stop('Weird!')
  }
  
  # Keep CN segments on chr of interest only
  # First, remove the 'chr' prefix in BTB segments if exists
  btb$chr[grepl('^chr',btb$chr)] = gsub('^chr','',btb$chr[grepl('^chr',btb$chr)])
  if(sum(btb$chr %in% tgtChrs) == 0){
    stop('No BTB segsments are found on targetted chromosomes. Please check to make sure that this is intentional')
  }
  
  btb = btb[btb$chr %in% tgtChrs,]
  btb$chr = ifelse(btb$chr == 'X',23,as.numeric(btb$chr))
  
  # convert to GRange object
  btb = GRanges(btb$chr,IRanges(btb$startpos,btb$endpos),
                chr=btb$chr,idx = btb$idx,posID=btb$posID,
                nMaj1_A=btb$nMaj1_A,nMin1_A=btb$nMin1_A, frac1_A=btb$frac1_A,
                nMaj2_A=btb$nMaj2_A,nMin2_A=btb$nMin2_A, frac2_A=btb$frac2_A)
  
  
  # Object containing final list of segments (without any subclonal segments?)
  clonalSegs = GRanges()
  # Object containing final list of sub_clonal segments only? 
  subClonalSegs = GRanges()
  
  for(chr in tgtChrs){
    #message(sprintf('Chromosome %s',chr))
    chr_segs = GRanges() # Object containing all segments on the chromosome
    data = btb[seqnames(btb) == chr]
    
    ## 1. Fill in the gaps - they should all be diploids
    gap_segs = GRanges(gaps(data),chr=chr,patNum=1, matNum=1,totCN=2,tumFrac = 0.5)
    # Fill in the end of the chromosome too (if needed)
    chromLen = max(chromInfo[chromInfo$chrom == chr,]$end)*1000
    if(max(end(data)) < chromLen){
      tmp = GRanges(chr,IRanges(max(end(data))+1,chromLen),chr=chr,patNum=1, matNum=1,totCN=2,tumFrac = 0.5)
      gap_segs = append(gap_segs,tmp)
    }
    # else if (max(end(data)) > chromLen){  # If max(end(data)) is longer than chromLen (just leave it for now...)
    #   if(start(data[end(data) == max(end(data))]) < chromLen){
    #     end(data[end(data) == max(end(data))]) = chromLen  
    #   }else{
    #     #remove that seg
    #     data = data[end(data) != max(end(data))]
    #   }
    #   
    # }
    
    chr_segs = append(chr_segs,gap_segs)
    
    
    ## 2. Extract all Clonal segments (i.e. segments which are present in ALL cells)
    clonalSegs = data[which(is.na(data$frac2_A))]
    if(length(clonalSegs)>0){
      mcols(clonalSegs) = data.frame(chr=chr,idx = clonalSegs$idx, patNum=clonalSegs$nMin1_A, matNum=clonalSegs$nMaj1_A,
                                     totCN=clonalSegs$nMin1_A+clonalSegs$nMaj1_A,tumFrac = clonalSegs$nMaj1_A/(clonalSegs$nMin1_A+clonalSegs$nMaj1_A))  
      if(any(clonalSegs$totCN == 0)){ # homozygous loss
        clonalSegs[clonalSegs$totCN == 0]$tumFrac = 0  
      }
      
    }
    
    
    ## 3. Extract all "non-clonal" segments (i.e. segments which are NOT present in ALL cells)
    sub_clonal = data[which((!is.na(data$frac2_A)))] 
    #sub_clonal = btb[which((!is.na(btb$frac2_A)) & (pmin(btb$frac1_A,btb$frac2_A) >= 0.1)),c('Idx','chr','startpos','endpos','nMaj1_A','nMin1_A','frac1_A','nMaj2_A','nMin2_A','frac2_A','posID','segLen')] 
    if(length(sub_clonal) > 0){
      message(sprintf('%s - chr %s: sub_clonal segments found!',PDID,chr))
      sub_clonal$segType = 'sub_clonal'  
      
      ### Sort out overlapping segments
      # There should be no overlaps between no_change and clonal / sub_clonal.
      # If there are overlaps between clonal and sub_clonal - remove the overlapping segments from clonal
      o = findOverlapPairs(clonalSegs,sub_clonal)
      # Remove overlapping segments from "clonalSegs" group (to be replaced by disjoin segments only)
      clonalSegs = clonalSegs[!clonalSegs$idx %in% o@first$idx]
      # # If there are any overlaping segments between "clonalSegs" and "sub_clonal"
      # if(length(o)>0){
      #   message(sprintf('%s - chr %s: OVERLAPPING sub_clonal segments found!',PDID,chr))
      #   non_clonalSegs_o = append(o@first,o@second)
      #   disjoin_ranges = disjoin(non_clonalSegs_o)
      #   disjoin_ranges = disjoin_ranges[!overlapsAny(disjoin_ranges,sub_clonal)]
      #   disjoin_ranges_non_clonalSegs = findOverlapPairs(disjoin_ranges,o@first)
      #   # Check that disjoin_ranges_non_clonalSegs@second should always be of type "clonalSegs"
      #   if(sum(disjoin_ranges_non_clonalSegs@second$type != 'clonalSegs') > 0){
      #     stop(sprintf('%s - chr %s: Something strange in disjoin happend...',PDID,chr))
      #   }
      #   # Each new clonalSegs segment should match to a unique segment from o3@first
      #   for(i in 1:length(disjoin_ranges_non_clonalSegs)){
      #     pair = disjoin_ranges_non_clonalSegs[i]
      #     start(pair@second) = start(pair@first)
      #     end(pair@second) = end(pair@first)
      #     clonalSegs = append(clonalSegs,pair@second)
      #   }
      # }
      
      # Restructure the sub_clonal segments
      
      for(i in 1:length(sub_clonal)){
        seg = sub_clonal[i]
        # Define maj vs sub clones
        if(seg$frac1_A > seg$frac2_A){
          maj_seg = GRanges(seg$chr,IRanges(start(seg),end(seg)),chr=seg$chr,clonalType='maj',
                            patNum=seg$nMin1_A, matNum=seg$nMaj1_A,totCN=seg$nMin1_A+seg$nMaj1_A,tumFrac = seg$nMaj1_A/(seg$nMaj1_A+seg$nMin1_A))
          sub_seg = GRanges(seg$chr,IRanges(start(seg),end(seg)),chr=seg$chr,clonalType='sub',
                            patNum=seg$nMin2_A, matNum=seg$nMaj2_A,totCN=seg$nMin2_A+seg$nMaj2_A,tumFrac = seg$nMaj2_A/(seg$nMaj2_A+seg$nMin2_A))
        }else{
          sub_seg = GRanges(seg$chr,IRanges(start(seg),end(seg)),chr=seg$chr,clonalType='sub',
                            patNum=seg$nMin1_A, matNum=seg$nMaj1_A,totCN=seg$nMin1_A+seg$nMaj1_A,tumFrac = seg$nMaj1_A/(seg$nMaj1_A+seg$nMin1_A))
          maj_seg = GRanges(seg$chr,IRanges(start(seg),end(seg)),chr=seg$chr,clonalType='maj',
                            patNum=seg$nMin2_A, matNum=seg$nMaj2_A,totCN=seg$nMin2_A+seg$nMaj2_A,tumFrac = seg$nMaj2_A/(seg$nMaj2_A+seg$nMin2_A))
        }
        subClonalSegs = append(subClonalSegs,sub_seg)
        subClonalSegs = append(subClonalSegs,maj_seg)
      }
      
    }
    
    
    chr_segs = append(chr_segs,clonalSegs)
    
    ### Merging segments of the same type and matNum:patNum configuration
    final_chr_segs = GRanges()
    chr_segs$tot2min = paste0(chr_segs$totCN,":",chr_segs$patNum)
    
    if(method == 'allelicRatio'){
      for(tumFrac in unique(chr_segs$tumFrac)){
        ### - MAYBE: If >= 95% of the chromosome has the same tumFrac (config) --> make the entire chromosome having the same config
        
        #print(config)
        chr_segs.sub = chr_segs[chr_segs$tumFrac == tumFrac]
        chr_segs.sub_ranges = GenomicRanges::reduce(chr_segs.sub)
        mcols(chr_segs.sub_ranges) = data.frame(chr=chr,tot2min=paste(unique(chr_segs.sub$tot2min),collapse = '_'),tumFrac = unique(chr_segs.sub$tumFrac))
        final_chr_segs = append(final_chr_segs,chr_segs.sub_ranges)
      }
    }else if(method == 'totalCN'){
      # Special scenarios for "PD46693" chr10 and 11
      #if((PDID == "PD46693" & chr %in% c(10,11)) |(PDID == "PD43255" & chr %in% c(18)) ){
      #  message(sprintf('%s chr%s: Replacing cn segments',PDID, chr))
      #  # Add artificial piece in
      #  final_chr_segs = GRanges(chr,IRanges(1,chromLen),chr=chr,tot2min='x:0',tumFrac = 1)
      #}
      for(totCN in unique(chr_segs$totCN)){
        ### - MAYBE: If >= 95% of the chromosome has the same tumFrac (config) --> make the entire chromosome having the same config
        
        #print(config)
        chr_segs.sub = chr_segs[chr_segs$totCN == totCN]
        chr_segs.sub_ranges = GenomicRanges::reduce(chr_segs.sub)
        mcols(chr_segs.sub_ranges) = data.frame(chr=chr,totCN=unique(chr_segs.sub$totCN),tumFrac = paste(unique(chr_segs.sub$tumFrac),collapse = '_'))
        final_chr_segs = append(final_chr_segs,chr_segs.sub_ranges)
      }
    }else if(method == 'config'){
      for(config in unique(chr_segs$tot2min)){
        chr_segs.sub = chr_segs[chr_segs$tot2min == config]
        chr_segs.sub_ranges = GenomicRanges::reduce(chr_segs.sub)
        totCN = as.numeric(strsplit(config,split=':')[[1]][1])
        patNum = as.numeric(strsplit(config,split=':')[[1]][2])
        matNum = totCN - patNum
        mcols(chr_segs.sub_ranges) = data.frame(chr=chr,totCN=totCN,matNum=matNum,patNum=patNum,tumFrac = paste(unique(chr_segs.sub$tumFrac),collapse = '_'))
        
        final_chr_segs = append(final_chr_segs,chr_segs.sub_ranges)  
      }
    }
    
    
    #final_chr_segs$tot2min = paste0(final_chr_segs$totCN,":",final_chr_segs$patNum)
    #final_chr_segs$posID = paste0(seqnames(final_chr_segs),':',start(final_chr_segs),'_',end(final_chr_segs))
    
    clonalSegs = append(clonalSegs,final_chr_segs)
  }
  
  
  
  ### Add subclone CNAs ###
  if(keepSubClonalSegs == T){
    clonalSegs$clonalType = 'none'
    # Remove segments of short length
    subClonalSegs = subClonalSegs[width(subClonalSegs) >= subCl.minSegLen]
    subClonalSegs$tot2min = paste0(subClonalSegs$totCN,":",subClonalSegs$patNum)
    finalSegs = append(clonalSegs,subClonalSegs)
  }else{
    finalSegs = clonalSegs
  }
  finalSegs$Start = start(finalSegs)
  finalSegs$Stop = end(finalSegs)
  #finalSegs$tot2min = paste0(finalSegs$totCN,":",finalSegs$patNum)
  finalSegs$posID = paste0(seqnames(finalSegs),':',start(finalSegs),'_',end(finalSegs))
  # Remove segments of short length
  finalSegs = finalSegs[width(finalSegs) >= minSegLen]
  
  
  
  
  if(longFormat == T){
    data = as.data.frame(mcols(finalSegs))
    out = pivot_longer(data,cols = c('Start','Stop'),names_to = 'posType',values_to = 'pos')
    out$idx = c(1:nrow(out))
    # Match max length for each chromosome
    for(chrom in unique(out$chr)){
      #current_maxChromLen = max(data[data$Chr == chrom,]$pos)
      maxChromLen = as.vector((chromInfo[(chromInfo$chrom == chrom) & (chromInfo$arm == 'q'),]$end)*1000)
      row = out[(out$chr == chrom) & (out$pos > maxChromLen),]$idx
      if(length(row) > 0){
        print(sprintf('Heey! %d %d',chrom,length(row)))
      }
    }
    
    # Get absolute genomic position
    out$abspos_kb = out$pos/1000 # if chromosome 1, abspos = pos
    for(r in 1:nrow(out)){
      chrom = out$chr[r]
      if (chrom > 1){
        out$abspos_kb[r] = out$abspos_kb[r] + (chromInfo[(chromInfo$chrom == (chrom-1)) & (chromInfo$arm == 'q'),]$abspos_kb)
      }
    }
    out$idx = c(1:nrow(out))
    out=out[order(out$chr,out$abspos_kb),]  
    if(removeBalancedSegs==T){
      out = out[out$tumFrac != 0.5,]
    }
  }else{
    out = finalSegs
    # Remove balanced segment
    if(removeBalancedSegs==T){
      out = out[out$tumFrac != 0.5]
    }
    out=out[order(out$chr,out$Start)]  
    out$idx = c(1:length(out))
    names(out) = seqnames(out)
  }
  
  return(out)
}










## version 2: Added automatic filter for inconsistent copy number segments
phaseSNPsFromCN_v2 = function (hSNPs, segs, refGenome, tBAM, outPath = NULL, FDR = 0.05, 
                               minPhasable = 2 * FDR, errRate = 0.01, useEM = TRUE, autoSegsFilter=TRUE,
                               alleleCounterParams = list(f = 3, F = 3852, m = 20, q = 35), nParallel = 1, plotMixtures = TRUE,
                               verbose = TRUE, ...) 
{
  if (is.null(names(segs))) 
    names(segs) = as.character(segs)
  if (any(duplicated(names(segs)))) 
    names(segs) = make.unique(names(segs))
  alleleCounterParams$x = FALSE
  alleleCounterParams$bams = tBAM
  alleleCounterParams$tgtLoci = hSNPs
  alleleCounterParams$refGenome = refGenome
  alleleCounterParams$outputs = outPath
  alleleCounterParams$nParallel = nParallel
  tCnts = do.call(alleleCounter, alleleCounterParams)[[1]]
  bases = c("A", "C", "G", "T")
  cnts = mcols(tCnts)
  tCnts$altCount = (as.matrix(cnts[, bases])[cbind(seq(nrow(cnts)), 
                                                   match(cnts$ALT, bases))])
  tCnts$refCount = (as.matrix(cnts[, bases])[cbind(seq(nrow(cnts)), 
                                                   match(cnts$REF, bases))])
  m = match(tCnts, hSNPs)
  hSNPs$altCountTum = hSNPs$refCountTum = hSNPs$totCountTum = hSNPs$altIsMum = hSNPs$matProb = hSNPs$passSanity = NA
  hSNPs$informative = FALSE
  hSNPs$altCountTum[m] = tCnts$altCount
  hSNPs$refCountTum[m] = tCnts$refCount
  hSNPs$totCountTum[m] = tCnts$Tot
  hSNPs$passSanity[m] = TRUE
  if (useEM) {
    inconsistentSegs = c()
    for (ii in seq_along(segs)) {
      if (verbose) 
        message(sprintf("Fitting mixture model for segment %s", 
                        names(segs[ii])))
      tmp = subsetByOverlaps(hSNPs, segs[ii])
      aC = tmp$altCountTum
      tC = tmp$totCountTum
      k = 2
      n = length(tC)
      rhos = runif(k)
      taus = rep(1/k, k)
      i = 0
      Q = -Inf
      while (TRUE) {
        ll = matrix(dbinom(rep(aC, k), rep(tC, k), rep(rhos, 
                                                       each = n), log = TRUE), nrow = n, ncol = k)
        states = t(t(ll) + log(taus))
        stateProbs = exp(states - apply(states, 1, max))
        stateProbs = stateProbs/rowSums(stateProbs)
        Qnew = sum(stateProbs * (rep(log(taus), each = n) + 
                                   ll))
        rhos = colSums(stateProbs * aC)/colSums(stateProbs * 
                                                  tC)
        taus = colMeans(stateProbs)
        dQ = (Qnew - Q)
        Q = Qnew
        if (verbose > 1) 
          message(sprintf("  %03d: rho=%g, tau=%g, dQ = %g", 
                          i, rhos[1], taus[1], dQ))
        i = i + 1
        if (abs(dQ) < 1e-06 | i == 5000) 
          break
      }
      o = order(rhos)
      rhos = rhos[o]
      taus = taus[o]
      stateProbs = stateProbs[, o]
      if (plotMixtures) {
        mixCols = c("#FF000080", "#0000FF80")
        breaks = seq(0, 1, length.out = 100)
        a = hist((aC/tC)[stateProbs[, 1] > 0.5], breaks = breaks, 
                 plot = FALSE)
        b = hist((aC/tC)[stateProbs[, 1] < 0.5], breaks = breaks, 
                 plot = FALSE)
        yMax = max(c(a$counts, b$counts))
        plot(a, border = FALSE, xlim = c(0, 1), col = mixCols[1], 
             xlab = "BAF", main = sprintf("seg %s, mixes %.02f (%.01f%%) & %.02f (%.01f%%)", 
                                          names(segs[ii]), rhos[1], 100 * taus[1], 
                                          rhos[2], 100 * taus[2]), ylim = c(0, yMax))
        plot(b, border = FALSE, xlim = c(0, 1), col = mixCols[2], 
             add = TRUE)
        abline(v = rhos, col = mixCols, lty = 2)
      }
      if (abs(sum(rhos) - 1) > 0.2 || max(taus) > 0.65) {
        if(autoSegsFilter){
          message(sprintf("Mixture solutions for CN segment %s inconsistent with CN change.  Found mixture with means %.02f (%.01f%%) and %.02f (%.01f%%)", 
                          names(segs[ii]), rhos[1], 100 * taus[1], rhos[2], 
                          100 * taus[2]))
          inconsistentSegs = c(inconsistentSegs,ii)  
        }else{
          stop(sprintf("Mixture solutions for CN segment %s inconsistent with CN change.  Found mixture with means %.02f (%.01f%%) and %.02f (%.01f%%)", 
                       names(segs[ii]), rhos[1], 100 * taus[1], rhos[2], 
                       100 * taus[2]))
        }
      }
      else if (verbose) {
        message("  Tumour DNA consistent with CN change at this segment.")
        mm = match(tmp, hSNPs)
        hSNPs$matProb[mm] = stateProbs[, 2]
      }
      
    }
    hSNPs$informative = !is.na(hSNPs$matProb) & (0.5 - abs(hSNPs$matProb - 0.5)) < FDR
    if(autoSegsFilter){
      message(sprintf("Removing %d segments where Tumour DNA is inconsistent with CN change provided.",length(inconsistentSegs)))
      # Remove inconsistent segments from "segs"
      segs = segs[!inconsistentSegs]  
    }
    
    
  }
  else {
    pVals = pbinom(pmax(tCnts$altCount, tCnts$refCount) - 
                     1, tCnts$Tot, 0.5, lower.tail = FALSE)
    qVals = p.adjust(pVals, method = "BH")
    hSNPs$matProb[m] = pVals
    hSNPs$informative[m] = qVals < FDR
  }
  if (errRate > 0) {
    pVals = pbinom(pmin(hSNPs$altCountTum, hSNPs$refCountTum) - 
                     1, hSNPs$totCountTum, errRate, lower.tail = FALSE)
    qVals = p.adjust(pVals, method = "BH")
    hSNPs$passSanity[which(qVals >= FDR)] = FALSE
  }
  hSNPs$altIsMum = ifelse(hSNPs$informative, hSNPs$altCountTum > 
                            hSNPs$refCountTum, NA)
  if (verbose) {
    message("###########\n# Summary #\n###########")
    message(sprintf("Of %s heterozygous SNPs:", prettyNum(length(hSNPs), 
                                                          big.mark = ",")))
  }
  o = findOverlaps(hSNPs, segs)
  for (i in seq(0, length(segs))) {
    qH = queryHits(o)[subjectHits(o) == i]
    lab = names(segs[i])
    if (i == 0) {
      qH = unique(queryHits(o))
      lab = "global"
    }
    nSNPs = length(qH)
    nPhased = sum(hSNPs$informative[qH], na.rm = TRUE)
    if (verbose) {
      message(sprintf("  %s (%.01f%%) fall in CN segment %s, of which:", 
                      prettyNum(nSNPs, big.mark = ","), nSNPs/length(hSNPs) * 
                        100, lab))
      message(sprintf("    %s (%.01f%%) could be phased, of which:", 
                      prettyNum(nPhased, big.mark = ","), nPhased/nSNPs * 
                        100))
      nSane = sum(hSNPs$passSanity[qH] & hSNPs$informative[qH], 
                  na.rm = TRUE)
      message(sprintf("      %s (%.01f%%) passed sanity checks, of which:", 
                      prettyNum(nSane, big.mark = ","), nSane/nPhased * 
                        100))
      nAltMum = sum(hSNPs$altIsMum[qH] & hSNPs$passSanity[qH] & 
                      hSNPs$informative[qH], na.rm = TRUE)
      message(sprintf("        %s (%.01f%%) have the maternal allele as ALT", 
                      prettyNum(nAltMum, big.mark = ","), nAltMum/nSane * 
                        100))
    }
    if (i != 0 && nPhased/nSNPs < minPhasable) {
      hSNPs$passSanity[qH] = FALSE
      hSNPs$altIsMum[qH] = NA
      warning(sprintf("Only %s of %s SNPs (%.02f%%) could be phased in copy number segment %s.  This segment may not contain a real CN change.", 
                      prettyNum(nPhased, big.mark = ","), prettyNum(nSNPs, 
                                                                    big.mark = ","), nPhased/nSNPs * 100, names(segs[i])))
    }
  }
  hSNPs$matAllele = ifelse(hSNPs$altIsMum, hSNPs$ALT, hSNPs$REF)
  hSNPs$patAllele = ifelse(hSNPs$altIsMum, hSNPs$REF, hSNPs$ALT)
  hSNPs$majorAllele = hSNPs$matAllele
  hSNPs$minorAllele = hSNPs$patAllele
  hSNPs@metadata$segs = segs
  return(hSNPs)
}











#' @param gCnts.segs filtered phCnts object
#' @param out_fp full file name and path to save the plot
#' @param nCell vector of number of cells per category
#' @param minRead only plot SNPs with this minimum coverage aggregated across all cells within the group
#' @param chrom chromosomes to plot
#' @param width width of the png
#' @param height height of the png

make_copyNumberPlot_cellClusters = function(gCnts.segs,out_fp,nCell,minRead = 10,chrom = paste0('chr',c(1:22)),doPlot=T,width=4000,height=1000){
  
  clusterCnts.segs = aggregateByLists(gCnts.segs, assays = c("altCount", "refCount"), gCnts.segs$clusterID, gCnts.segs$regionID)
  clusterCnts.segs$pos = as.numeric(gsub('.*:|_.*$','',clusterCnts.segs$regionID))
  clusterCnts.segs$chr = gsub(':.*$','',clusterCnts.segs$regionID)
  clusterCnts.segs$altFreq = clusterCnts.segs$altCount / (clusterCnts.segs$altCount + clusterCnts.segs$refCount)
  clusterCnts.segs$totCount = clusterCnts.segs$altCount + clusterCnts.segs$refCount
  clusterCnts.segs$cov = ifelse(clusterCnts.segs$totCount < 5,'1-5',
                                ifelse(clusterCnts.segs$totCount < 10,'5-10',
                                       ifelse(clusterCnts.segs$totCount < 20,'10-20','>=20')))
  clusterCnts.segs$cov = factor(clusterCnts.segs$cov,c('1-5','5-10','10-20','>=20'))
  
  colnames(clusterCnts.segs)[colnames(clusterCnts.segs) == 'cellID'] = 'clusterID'
  
  # Count number of cells
  clusterCnts.segs$clusterID_nCell = as.vector(nCell[match(clusterCnts.segs$clusterID,names(nCell))])
  clusterCnts.segs$chr = factor(clusterCnts.segs$chr,paste0('chr',c(1:22,'X')))
  clusterCnts.segs$clusterID = paste0(clusterCnts.segs$clusterID,' (n=',clusterCnts.segs$clusterID_nCell,')')
  
  
  df = clusterCnts.segs[clusterCnts.segs$totCount > minRead & clusterCnts.segs$chr %in% chrom,]
  
  ## check if there is phasing information
  if('altIsMum' %in% colnames(mcols(gCnts_wUninf))){
    df$phasingAssign = gCnts_wUninf$altIsMum[match(df$regionID,names(gCnts_wUninf))]
    table(df$phasingAssign)
    df$phasingAssign = ifelse(is.na(df$phasingAssign), 'uninformative',
                              ifelse(df$phasingAssign == T,'major_allele','minor_allele'))
    
    
    ccs = c('major_allele' = col25[2],
            'minor_allele' = col25[1],
            'uninformative' = grey(0.8))  
    
    p = ggplot(df,aes(pos/1e6,altFreq))+
      geom_point(aes(col=phasingAssign),size=0.3,alpha=0.6)+
      facet_grid(clusterID  ~ chr,scales = 'free_x')+
      geom_hline(yintercept = 0.5,lty=2,lwd=0.5,col='black')+
      scale_color_manual(values = ccs)+
      scale_y_continuous(breaks = c(0.0,0.5,1.0),labels = c(0,'',1))+
      theme_classic(base_size = 14)+
      theme(panel.border = element_rect(fill=F),
            axis.line = element_blank(),
            axis.text = element_text(size=11),strip.text.y.right = element_text(size=9))+
      xlab('Genomic position') + ylab('Aggregated Alt allele frequency')
    
    
  }else{
    
    p = ggplot(df,aes(pos/1e6,altFreq))+
      geom_point(size=0.3,alpha=0.6)+
      facet_grid(clusterID  ~ chr,scales = 'free_x')+
      geom_hline(yintercept = 0.5,lty=2,lwd=0.5,col='black')+
      scale_color_manual(values = grey(0.6))+
      scale_y_continuous(breaks = c(0.0,0.5,1.0),labels = c(0,'',1))+
      theme_classic(base_size = 14)+
      theme(panel.border = element_rect(fill=F),
            axis.line = element_blank(),
            axis.text = element_text(size=11),strip.text.y.right = element_text(size=9))+
      xlab('Genomic position') + ylab('Aggregated Alt allele frequency')
  }
  
  
  x_max = max(df$pos)/1e6
  
  
  
  print(p)
  
  
  if(doPlot){
    png(out_fp,width = width,height = height)
    print(p)
    dev.off()  
  }
  
  return(df)
}



