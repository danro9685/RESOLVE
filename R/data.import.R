#' Create Single Base Substitutions (SBS) counts matrix from input data for a provided reference genome.
#'
#' @examples
#' library('BSgenome.Hsapiens.1000genomes.hs37d5')
#' data(ssm560_reduced)
#' res <- getSBSCounts(data = ssm560_reduced, reference = BSgenome.Hsapiens.1000genomes.hs37d5)
#'
#' @title getSBSCounts
#' @param data A data.frame with variants having 6 columns: sample name, chromosome, start position, end position, ref, alt.
#' @param reference A BSgenome object with the reference genome to be used to retrieve flanking bases.
#' @return A matrix with Single Base Substitutions (SBS) counts per patient.
#' @export getSBSCounts
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import IRanges
#' @importFrom data.table data.table dcast .N
#' @importFrom Biostrings DNAStringSet complement reverseComplement subseq
#' @importFrom BSgenome getSeq
#'
getSBSCounts <- function(data, reference = NULL) {

    # check that reference is a BSgenome object
    if (is.null(reference) | (!is(reference, "BSgenome"))) {
        stop("The reference genome provided as input needs to be a BSgenome object.")
    }

    # preprocessing input data
    data <- as.data.frame(data)
    colnames(data) <- c("sample", "chrom", "start", "end", "ref", "alt")

    # consider only single nucleotide variants involving (A,C,G,T) bases
    data <- data[which(data[, "start"] == data[, "end"]), , drop = FALSE]
    data <- data[which(as.character(data[, "ref"]) %in% c("A", "C", "G", "T")), ,
        drop = FALSE]
    data <- data[which(as.character(data[, "alt"]) %in% c("A", "C", "G", "T")), ,
        drop = FALSE]
    data <- data[, c("sample", "chrom", "start", "ref", "alt"), drop = FALSE]
    colnames(data) <- c("sample", "chrom", "pos", "ref", "alt")
    data <- unique(data)
    data <- data[order(data[, "sample"], data[, "chrom"], data[, "pos"]), , drop = FALSE]

    # convert data to GRanges
    data <- GRanges(data$chrom, IRanges(start = (data$pos - 1), width = 3), ref = DNAStringSet(data$ref),
        alt = DNAStringSet(data$alt), sample = data$sample)

    # check that all chromosomes match reference
    if (length(setdiff(seqnames(data), seqnames(reference))) > 0) {
        warning("Check chromosome names, not all match reference genome.")
    }

    # find context for each mutation
    data$context <- getSeq(reference, data)

    # check for any mismatch with BSgenome context
    if (any(subseq(data$context, 2, 2) != data$ref)) {
        warning("Check reference bases, not all match context.")
    }

    # get complements and reverse complements
    data$cref <- complement(data$ref)
    data$calt <- complement(data$alt)
    data$rccontext <- reverseComplement(data$context)

    # identify trinucleotides motif
    data$cat <- ifelse(as.character(data$ref) %in% c("C", "T"), paste0(subseq(data$context,
        1, 1), "[", data$ref, ">", data$alt, "]", subseq(data$context, 3, 3)), paste0(subseq(data$rccontext,
        1, 1), "[", data$cref, ">", data$calt, "]", subseq(data$rccontext, 3, 3)))

    # create 96 trinucleotides mutation categories
    categories_context <- NULL
    categories_alt <- rep(c(rep("C>A", 4), rep("C>G", 4), rep("C>T", 4), rep("T>A",
        4), rep("T>C", 4), rep("T>G", 4)), 4)
    categories_cat <- NULL
    cont <- 0
    for (i in c("A", "C", "G", "T")) {
        for (j in seq_len(6)) {
            for (k in c("A", "C", "G", "T")) {
                cont <- cont + 1
                categories_context <- c(categories_context, paste0(k, ":", i))
                categories_cat <- c(categories_cat, paste0(k, "[", categories_alt[cont],
                  "]", i))
            }
        }
    }
    mutation_categories <- data.table(context = categories_context, alt = categories_alt,
        cat = categories_cat)

    # count number of mutations per sample for each category
    data <- merge(mutation_categories[, list(cat)], data.table(sample = data$sample,
        cat = data$cat)[, .N, by = list(sample, cat)], by = "cat", all = TRUE)
    data <- dcast(data, sample ~ cat, value.var = "N")
    data <- data[!is.na(sample), drop = FALSE]
    data[is.na(data)] <- 0

    # make trinucleotides counts matrix
    samples_names <- data$sample
    data <- as.matrix(data[, 2:ncol(data), drop = FALSE])
    rownames(data) <- samples_names
    data <- data[sort(rownames(data)), , drop = FALSE]
    data <- data[, sort(colnames(data)), drop = FALSE]
    trinucleotides_counts <- array(0, c(nrow(data), 96))
    rownames(trinucleotides_counts) <- rownames(data)
    colnames(trinucleotides_counts) <- sort(as.character(mutation_categories$cat))
    rows_contexts <- rownames(data)
    cols_contexts <- colnames(trinucleotides_counts)[which(colnames(trinucleotides_counts) %in%
        colnames(data))]
    trinucleotides_counts[rows_contexts, cols_contexts] <- data[rows_contexts, cols_contexts]

    # return trinucleotides counts matrix
    return(trinucleotides_counts)

}

#' Create Multi-Nucleotide Variants (MNVs) counts matrix from input data.
#'
#' @examples
#' data(ssm560_reduced)
#' res <- getMNVCounts(data = ssm560_reduced)
#'
#' @title getMNVCounts
#' @param data A data.frame with variants having 6 columns: sample name, chromosome, start position, end position, ref, alt.
#' @param predefined_dbs_mbs Boolean. As defined by the function get_mut_type from the package MutationalPatterns, it specifies whether 
#' dbs and mbs mutations have been predefined in the input data. This function by default assumes that dbs and mbs mutations are present 
#' in the vcf as snvs, which are positioned next to each other. If your dbs/mbs mutations are called separately, you should set this 
#' argument to TRUE.
#' @return A matrix with Multi-Nucleotide Variants (MNVs) counts per patient.
#' @export getMNVCounts
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import IRanges
#' @importFrom MutationalPatterns get_mut_type get_dbs_context count_dbs_contexts
#'
getMNVCounts <- function( data, predefined_dbs_mbs = FALSE ) {

    # preprocessing input data
    data <- as.data.frame(data)
    colnames(data) <- c("sample", "chrom", "start", "end", "ref", "alt")
    data <- data[, c("sample", "chrom", "start", "ref", "alt"), drop = FALSE]
    colnames(data) <- c("sample", "chrom", "pos", "ref", "alt")
    data <- unique(data)
    data <- data[order(data[, "sample"], data[, "chrom"], data[, "pos"]), , drop = FALSE]

    # convert data to GRanges
    data <- GRanges(data$chrom, IRanges(start = data$pos, width = nchar(data$ref)), 
        ref = data$ref, alt = data$alt, sample = data$sample)
    data <- get_mut_type(data, type = "dbs", predefined_dbs_mbs = predefined_dbs_mbs)
    data <- get_dbs_context(data)

    # build counts matrix
    dbs_counts <- NULL
    samp_names <- sort(unique(data$sample))
    for(i in samp_names) {
        res <- t(count_dbs_contexts(data[which(data$sample==i), , drop = FALSE]))
        dbs_counts <- rbind(dbs_counts, res)
    }
    rownames(dbs_counts) <- samp_names

    # create 78 doublet base substitutions categories and their complements
    mutation_categories <- c("AC>CA", "AC>CG", "AC>CT", "AC>GA", "AC>GG", "AC>GT",
        "AC>TA", "AC>TG", "AC>TT")  # AC
    mutation_categories <- c(mutation_categories, "AT>CA", "AT>CC", "AT>CG", "AT>GA",
        "AT>GC", "AT>TA")  # AT
    mutation_categories <- c(mutation_categories, "CC>AA", "CC>AG", "CC>AT", "CC>GA",
        "CC>GG", "CC>GT", "CC>TA", "CC>TG", "CC>TT")  # CC
    mutation_categories <- c(mutation_categories, "CG>AT", "CG>GC", "CG>GT", "CG>TA",
        "CG>TC", "CG>TT")  # CG
    mutation_categories <- c(mutation_categories, "CT>AA", "CT>AC", "CT>AG", "CT>GA",
        "CT>GC", "CT>GG", "CT>TA", "CT>TC", "CT>TG")  # CT
    mutation_categories <- c(mutation_categories, "GC>AA", "GC>AG", "GC>AT", "GC>CA",
        "GC>CG", "GC>TA")  # GC
    mutation_categories <- c(mutation_categories, "TA>AT", "TA>CG", "TA>CT", "TA>GC",
        "TA>GG", "TA>GT")  # TA
    mutation_categories <- c(mutation_categories, "TC>AA", "TC>AG", "TC>AT", "TC>CA",
        "TC>CG", "TC>CT", "TC>GA", "TC>GG", "TC>GT")  # TC
    mutation_categories <- c(mutation_categories, "TG>AA", "TG>AC", "TG>AT", "TG>CA",
        "TG>CC", "TG>CT", "TG>GA", "TG>GC", "TG>GT")  # TG
    mutation_categories <- c(mutation_categories, "TT>AA", "TT>AC", "TT>AG", "TT>CA",
        "TT>CC", "TT>CG", "TT>GA", "TT>GC", "TT>GG")  # TT
    mutation_categories_complements <- c("GT>TG", "GT>CG", "GT>AG", "GT>TC", "GT>CC",
        "GT>AC", "GT>TA", "GT>CA", "GT>AA")  # GT
    mutation_categories_complements <- c(mutation_categories_complements, "AT>TG",
        "AT>GG", "AT>CG", "AT>TC", "AT>GC", "AT>TA")  # AT
    mutation_categories_complements <- c(mutation_categories_complements, "GG>TT",
        "GG>CT", "GG>AT", "GG>TC", "GG>CC", "GG>AC", "GG>TA", "GG>CA", "GG>AA")  # GG
    mutation_categories_complements <- c(mutation_categories_complements, "CG>AT",
        "CG>GC", "CG>AC", "CG>TA", "CG>GA", "CG>AA")  # CG
    mutation_categories_complements <- c(mutation_categories_complements, "AG>TT",
        "AG>GT", "AG>CT", "AG>TC", "AG>GC", "AG>CC", "AG>TA", "AG>GA", "AG>CA")  # AG
    mutation_categories_complements <- c(mutation_categories_complements, "GC>TT",
        "GC>CT", "GC>AT", "GC>TG", "GC>CG", "GC>TA")  # GC
    mutation_categories_complements <- c(mutation_categories_complements, "TA>AT",
        "TA>CG", "TA>AG", "TA>GC", "TA>CC", "TA>AC")  # TA
    mutation_categories_complements <- c(mutation_categories_complements, "GA>TT",
        "GA>CT", "GA>AT", "GA>TG", "GA>CG", "GA>AG", "GA>TC", "GA>CC", "GA>AC")  # GA
    mutation_categories_complements <- c(mutation_categories_complements, "CA>TT",
        "CA>GT", "CA>AT", "CA>TG", "CA>GG", "CA>AG", "CA>TC", "CA>GC", "CA>AC")  # CA
    mutation_categories_complements <- c(mutation_categories_complements, "AA>TT",
        "AA>GT", "AA>CT", "AA>TG", "AA>GG", "AA>CG", "AA>TC", "AA>GC", "AA>CC")  # AA
    colnames(dbs_counts) = gsub("_",">",colnames(dbs_counts))

    # make the multi nucleotide counts matrix
    multi_nucleotides_counts = dbs_counts[, mutation_categories ,drop = FALSE]

    # return the multi nucleotides counts matrix
    return(multi_nucleotides_counts)

}

#' Create Small Insertions and Deletions (IDs) counts matrix from input data.
#'
#' @examples
#' library('BSgenome.Hsapiens.1000genomes.hs37d5')
#' data(id_example_reduced)
#' res <- getIDCounts(data = id_example_reduced, reference = BSgenome.Hsapiens.1000genomes.hs37d5)
#'
#' @title getIDCounts
#' @param data A data.frame with variants having 6 columns: sample name, chromosome, start position, end position, ref, alt.
#' @param reference A BSgenome object with the reference genome to be used.
#' @return A matrix with Small Insertions and Deletions (IDs) counts per patient.
#' @export getIDCounts
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import IRanges
#' @importFrom MutationalPatterns get_mut_type get_indel_context count_indel_contexts
#' @importFrom S4Vectors metadata
#'
getIDCounts <- function(data, reference = NULL) {

    # check that reference is a BSgenome object
    if (is.null(reference) | (!is(reference, "BSgenome"))) {
        stop("The reference genome provided as input needs to be a BSgenome object.")
    }

    # preprocessing input data
    data <- as.data.frame(data)
    colnames(data) <- c("sample", "chrom", "start", "end", "ref", "alt")
    data <- data[, c("sample", "chrom", "start", "ref", "alt"), drop = FALSE]
    colnames(data) <- c("sample", "chrom", "pos", "ref", "alt")
    data <- unique(data)
    data <- data[order(data[, "sample"], data[, "chrom"], data[, "pos"]), , drop = FALSE]

    # convert data to GRanges
    data <- GRanges(data$chrom, IRanges(start = data$pos, width = nchar(data$ref)), 
        ref = data$ref, alt = data$alt, sample = data$sample)
    data <- get_mut_type(data, type = "indel", predefined_dbs_mbs = FALSE)
    genome(data) <- metadata(reference)$genome
    data <- get_indel_context(data, reference)

    # build counts matrix
    id_counts <- NULL
    samp_names <- sort(unique(data$sample))
    for(i in samp_names) {
        res <- t(count_indel_contexts(data[which(data$sample==i), , drop = FALSE]))
        id_counts <- rbind(id_counts, res)
    }
    rownames(id_counts) <- samp_names

    # rename mutation categories
    mutation_categories <- c("1:Del:C:0" , "1:Del:C:1" , "1:Del:C:2" , "1:Del:C:3" , "1:Del:C:4" , "1:Del:C:5" , 
        "1:Del:T:0" , "1:Del:T:1" , "1:Del:T:2" , "1:Del:T:3" , "1:Del:T:4" , "1:Del:T:5" , "1:Ins:C:0" , "1:Ins:C:1" , 
        "1:Ins:C:2" , "1:Ins:C:3" , "1:Ins:C:4" , "1:Ins:C:5" , "1:Ins:T:0" , "1:Ins:T:1" , "1:Ins:T:2" , "1:Ins:T:3" , 
        "1:Ins:T:4" , "1:Ins:T:5" , "2:Del:R:0" , "2:Del:R:1" , "2:Del:R:2" , "2:Del:R:3" , "2:Del:R:4" , "2:Del:R:5" , 
        "3:Del:R:0" , "3:Del:R:1" , "3:Del:R:2" , "3:Del:R:3" , "3:Del:R:4" , "3:Del:R:5" , "4:Del:R:0" , "4:Del:R:1" , 
        "4:Del:R:2" , "4:Del:R:3" , "4:Del:R:4" , "4:Del:R:5" , "5:Del:R:0" , "5:Del:R:1" , "5:Del:R:2" , "5:Del:R:3" , 
        "5:Del:R:4" , "5:Del:R:5" , "2:Ins:R:0" , "2:Ins:R:1" , "2:Ins:R:2" , "2:Ins:R:3" , "2:Ins:R:4" , "2:Ins:R:5" , 
        "3:Ins:R:0" , "3:Ins:R:1" , "3:Ins:R:2" , "3:Ins:R:3" , "3:Ins:R:4" , "3:Ins:R:5" , "4:Ins:R:0" , "4:Ins:R:1" , 
        "4:Ins:R:2" , "4:Ins:R:3" , "4:Ins:R:4" , "4:Ins:R:5" , "5:Ins:R:0" , "5:Ins:R:1" , "5:Ins:R:2" , "5:Ins:R:3" , 
        "5:Ins:R:4" , "5:Ins:R:5" , "2:Del:M:1" , "3:Del:M:1" , "3:Del:M:2" , "4:Del:M:1" , "4:Del:M:2" , "4:Del:M:3" , 
        "5:Del:M:1" , "5:Del:M:2" , "5:Del:M:3" , "5:Del:M:4" , "5:Del:M:5")
    colnames(id_counts) <- mutation_categories

    # return indel counts matrix
    return(id_counts)

}

#' Create Copy Numbers (CNs) counts matrix from input data.
#' This function has been derived from: https://github.com/UCL-Research-Department-of-Pathology/panConusig/blob/main/R/setup_CNsigs.R
#'
#' @examples
#' data(cn_example_reduced)
#' res <- getCNCounts(data = cn_example_reduced)
#'
#' @title getCNCounts
#' @param data A data.frame with copy number data having 6 columns: sample name, chromosome, start position, end position, major CN, minor CN.
#' @return A matrix with Copy Numbers (CNs) counts per patient.
#' @export getCNCounts
#'
getCNCounts <- function( data ) {

    # setup all possible CN classes
    CNclasses <- c("1","2","3-4","5-8","9+") # different total CN states
    lengths <- c("0-100kb","100kb-1Mb","1Mb-10Mb","10Mb-40Mb",">40Mb") # different lengths
    lohs <- c("LOH","het") # loh statuses
    allNames <- sapply(lohs,FUN=function(x) sapply(CNclasses,FUN=function(y) sapply(lengths,FUN=function(z) paste0(y,":",x,":",z))))
    allNames <- allNames[-grep("1[:]het",allNames)] # remove heterozygous deletions (impossible)
    homdelclasses <- "0" # add homozygous deletions
    homdelLengths <- c("0-100kb","100kb-1Mb",">1Mb")
    allNames <- c(sapply(homdelclasses,FUN=function(x) sapply(homdelLengths,FUN=function(z) paste0(x,":homdel:",z))),allNames)

    # set colnames for data
    colnames(data) <- c("sample","chrom","start","end","major","minor")

    # compute CN data for each patients
    sample_names <- sort(unique(data$sample))
    CN_data <- sapply(sample_names, FUN = function(x) {
        patient_pos <- which(data$sample==x)
        chrom <- as.character(data$chrom[patient_pos])
        start <- as.numeric(data$start[patient_pos])
        end <- as.numeric(data$end[patient_pos])
        major <- as.numeric(data$major[patient_pos])
        minor <- as.numeric(data$minor[patient_pos])
        total_cn <- (major+minor)
        lengths <- (end-start)/1000000
        LOH <- pmin(major,minor)
        LOHstatus <- ifelse(LOH==0,"LOH","het")
        LOHstatus[which(total_cn==0)] <- "homdel"
        varCN <- cut(total_cn, 
                    breaks = c(-0.5,0.5,1.5,2.5,4.5,8.5,Inf), 
                    labels = c("0","1","2","3-4","5-8","9+"))
        varLength <- rep(NA, length = length(varCN))
        hdIndex <- which(LOHstatus=="homdel")
        if(length(hdIndex)>0) {
            varLength[hdIndex] <- paste0(cut(lengths[hdIndex],breaks=c(-0.01,0.1,1,Inf)))
            varLength[-hdIndex] <- paste0(cut(lengths[-hdIndex],breaks=c(-0.01,0.1,1,10,40,Inf)))
        } else {
            varLength <- paste0(cut(lengths,breaks=c(-0.01,0.1,1,10,40,Inf)))
        }
        renameVarLengths <- c("(-0.01,0.1]"="0-100kb","(0.1,1]"="100kb-1Mb","(1,10]"="1Mb-10Mb","(10,40]"="10Mb-40Mb","(40,Inf]"=">40Mb","(1,Inf]"=">1Mb")
        varLength <- renameVarLengths[paste0(varLength)]
        sepVars <- paste(varCN, LOHstatus, varLength, sep=":")
        variables <- table(sepVars)
    })
    names(CN_data) <- sample_names

    # build CN counts matrix
    cn_counts <- matrix(0, nrow = length(CN_data), ncol = length(allNames))
    rownames(cn_counts) <- names(CN_data)
    colnames(cn_counts) <- allNames
    for(i in rownames(cn_counts)) {
        cn_counts[i,names(CN_data[[i]])] <- CN_data[[i]]
    }

    # return copy numbers counts matrix
    return(cn_counts)

}
