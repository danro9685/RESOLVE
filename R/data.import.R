#' Create Single Base Substitutions (SBS) counts matrix from input data for a provided reference genome.
#'
#' @examples
#' library("BSgenome.Hsapiens.1000genomes.hs37d5")
#' data(ssm560_reduced)
#' res = getSBSCounts(data = ssm560_reduced, reference = BSgenome.Hsapiens.1000genomes.hs37d5)
#'
#' @title getSBSCounts
#' @param data a data.frame with variants having 6 columns: sample name, chromosome, start position, end position, ref, alt.
#' @param reference a BSgenome object with the reference genome to be used to retrieve flanking bases.
#' @return A matrix with Single Base Substitutions (SBS) counts per patient.
#' @export getSBSCounts
#' @import IRanges
#' @import GenomeInfoDb
#' @import BSgenome.Hsapiens.1000genomes.hs37d5
#' @importFrom data.table data.table dcast .N
#' @importFrom Biostrings DNAStringSet complement reverseComplement subseq
#' @importFrom BSgenome getSeq
#' @importFrom GenomicRanges GRanges seqnames
#'
"getSBSCounts" <- function( data, reference = NULL ) {

    # check that reference is a BSgenome object
    if(is.null(reference)|(!is(reference,"BSgenome"))) {
        stop("The reference genome provided as input needs to be a BSgenome object.")
    }

    # preprocessing input data
    data <- as.data.frame(data)
    colnames(data) <- c("sample","chrom","start","end","ref","alt")

    # consider only single nucleotide variants involving (A,C,G,T) bases
    data <- data[which(data[,"start"]==data[,"end"]),,drop=FALSE]
    data <- data[which(as.character(data[,"ref"])%in%c("A","C","G","T")),,drop=FALSE]
    data <- data[which(as.character(data[,"alt"])%in%c("A","C","G","T")),,drop=FALSE]
    data <- data[,c("sample","chrom","start","ref","alt"),drop=FALSE]
    colnames(data) <- c("sample","chrom","pos","ref","alt")
    data <- unique(data)
    data <- data[order(data[,"sample"],data[,"chrom"],data[,"pos"]),,drop=FALSE]

    # convert data to GRanges
    data <- GRanges(data$chrom,IRanges(start=(data$pos-1),width=3),ref=DNAStringSet(data$ref),alt=DNAStringSet(data$alt),sample=data$sample)

    # check that all chromosomes match reference
    if(length(setdiff(seqnames(data),seqnames(reference)))>0) {
        warning("Check chromosome names, not all match reference genome.")
    }

    # find context for each mutation
    data$context <- getSeq(reference,data)

    # check for any mismatch with BSgenome context
    if(any(subseq(data$context,2,2)!=data$ref)) {
        warning("Check reference bases, not all match context.")
    }

    # get complements and reverse complements
    data$cref <- complement(data$ref)
    data$calt <- complement(data$alt)
    data$rccontext <- reverseComplement(data$context)

    # identify trinucleotides motif
    data$cat <- ifelse(as.character(data$ref)%in%c("C","T"),paste0(subseq(data$context,1,1),"[",data$ref,">",data$alt,"]",subseq(data$context,3,3)), 
                       paste0(subseq(data$rccontext,1,1),"[",data$cref,">",data$calt,"]",subseq(data$rccontext,3,3)))

    # create 96 trinucleotides mutation categories
    categories_context <- NULL
    categories_alt <- rep(c(rep("C>A",4),rep("C>G",4),rep("C>T",4),rep("T>A",4),rep("T>C",4),rep("T>G",4)),4)
    categories_cat <- NULL
    cont <- 0
    for(i in c("A","C","G","T")) {
        for(j in 1:6) {
            for(k in c("A","C","G","T")) {
                cont <- cont + 1
                categories_context <- c(categories_context,paste0(k,":",i))
                categories_cat <- c(categories_cat,paste0(k,"[",categories_alt[cont],"]",i))
            }
        }
    }
    mutation_categories <- data.table(context=categories_context,alt=categories_alt,cat=categories_cat)
    
    # count number of mutations per sample for each category
    data <- merge(mutation_categories[,list(cat)],data.table(sample=data$sample,cat=data$cat)[,.N,by=list(sample,cat)],by="cat",all=TRUE)
    data <- dcast(data,sample~cat,value.var="N")
    data <- data[!is.na(sample),drop=FALSE]
    data[is.na(data)] <- 0

    # make trinucleotides counts matrix
    samples_names <- data$sample
    data <- as.matrix(data[,2:ncol(data),drop=FALSE])
    rownames(data) <- samples_names
    data <- data[sort(rownames(data)),,drop=FALSE]
    data <- data[,sort(colnames(data)),drop=FALSE]
    trinucleotides_counts <- array(0,c(nrow(data),96))
    rownames(trinucleotides_counts) <- rownames(data)
    colnames(trinucleotides_counts) <- sort(as.character(mutation_categories$cat))
    rows_contexts <- rownames(data)
    cols_contexts <- colnames(trinucleotides_counts)[which(colnames(trinucleotides_counts)%in%colnames(data))]
    trinucleotides_counts[rows_contexts,cols_contexts] <- data[rows_contexts,cols_contexts]

    # return trinucleotides counts matrix
    return(trinucleotides_counts)

}

#' Create Multi-Nucleotide Variants (MNVs) counts matrix from input data.
#'
#' @examples
#' data(ssm560_reduced)
#' res = getMNVCounts(data = ssm560_reduced)
#'
#' @title getMNVCounts
#' @param data a data.frame with variants having 6 columns: sample name, chromosome, start position, end position, ref, alt.
#' @return A matrix with Multi-Nucleotide Variants (MNVs) counts per patient.
#' @export getMNVCounts
#' @importFrom data.table data.table dcast .N
#'
"getMNVCounts" <- function( data ) {

    # preprocessing input data
    data <- as.data.frame(data)
    colnames(data) <- c("sample","chrom","start","end","ref","alt")

    # consider only single nucleotide variants involving (A,C,G,T) bases
    data <- data[which(data[,"start"]==data[,"end"]),,drop=FALSE]
    data <- data[which(as.character(data[,"ref"])%in%c("A","C","G","T")),,drop=FALSE]
    data <- data[which(as.character(data[,"alt"])%in%c("A","C","G","T")),,drop=FALSE]
    data <- data[,c("sample","chrom","start","ref","alt"),drop=FALSE]
    colnames(data) <- c("sample","chrom","pos","ref","alt")
    data <- unique(data)
    data <- data[order(data[,"sample"],data[,"chrom"],data[,"pos"]),,drop=FALSE]

    # consider Multi-Nucleotide Variants (MNVs)
    cond1 <- c(data$sample[1:(nrow(data)-1)]==data$sample[2:nrow(data)],FALSE)
    cond2 <- c(data$chrom[1:(nrow(data)-1)]==data$chrom[2:nrow(data)],FALSE)
    cond3 <- c(data$pos[1:(nrow(data)-1)]==(data$pos[2:nrow(data)]-1),FALSE)
    data <- cbind(data[which(cond1&cond2&cond3),],data[(which(cond1&cond2&cond3)+1),c("ref","alt")])
    data[,4] <- paste0(data[,4],data[,6])
    data[,5] <- paste0(data[,5],data[,7])
    data <- data[,c(1:5),drop=FALSE]

    # create 78 doublet base substitutions categories and their complements
    mutation_categories <- c("AC>CA","AC>CG","AC>CT","AC>GA","AC>GG","AC>GT","AC>TA","AC>TG","AC>TT") # AC
    mutation_categories <- c(mutation_categories,"AT>CA","AT>CC","AT>CG","AT>GA","AT>GC","AT>TA") # AT
    mutation_categories <- c(mutation_categories,"CC>AA","CC>AG","CC>AT","CC>GA","CC>GG","CC>GT","CC>TA","CC>TG","CC>TT") # CC
    mutation_categories <- c(mutation_categories,"CG>AT","CG>GC","CG>GT","CG>TA","CG>TC","CG>TT") # CG
    mutation_categories <- c(mutation_categories,"CT>AA","CT>AC","CT>AG","CT>GA","CT>GC","CT>GG","CT>TA","CT>TC","CT>TG") # CT
    mutation_categories <- c(mutation_categories,"GC>AA","GC>AG","GC>AT","GC>CA","GC>CG","GC>TA") # GC
    mutation_categories <- c(mutation_categories,"TA>AT","TA>CG","TA>CT","TA>GC","TA>GG","TA>GT") # TA
    mutation_categories <- c(mutation_categories,"TC>AA","TC>AG","TC>AT","TC>CA","TC>CG","TC>CT","TC>GA","TC>GG","TC>GT") # TC
    mutation_categories <- c(mutation_categories,"TG>AA","TG>AC","TG>AT","TG>CA","TG>CC","TG>CT","TG>GA","TG>GC","TG>GT") # TG
    mutation_categories <- c(mutation_categories,"TT>AA","TT>AC","TT>AG","TT>CA","TT>CC","TT>CG","TT>GA","TT>GC","TT>GG") # TT
    mutation_categories_complements <- c("GT>TG","GT>CG","GT>AG","GT>TC","GT>CC","GT>AC","GT>TA","GT>CA","GT>AA") # GT
    mutation_categories_complements <- c(mutation_categories_complements,"AT>TG","AT>GG","AT>CG","AT>TC","AT>GC","AT>TA") # AT
    mutation_categories_complements <- c(mutation_categories_complements,"GG>TT","GG>CT","GG>AT","GG>TC","GG>CC","GG>AC","GG>TA","GG>CA","GG>AA") # GG
    mutation_categories_complements <- c(mutation_categories_complements,"CG>AT","CG>GC","CG>AC","CG>TA","CG>GA","CG>AA") # CG
    mutation_categories_complements <- c(mutation_categories_complements,"AG>TT","AG>GT","AG>CT","AG>TC","AG>GC","AG>CC","AG>TA","AG>GA","AG>CA") # AG
    mutation_categories_complements <- c(mutation_categories_complements,"GC>TT","GC>CT","GC>AT","GC>TG","GC>CG","GC>TA") # GC
    mutation_categories_complements <- c(mutation_categories_complements,"TA>AT","TA>CG","TA>AG","TA>GC","TA>CC","TA>AC") # TA
    mutation_categories_complements <- c(mutation_categories_complements,"GA>TT","GA>CT","GA>AT","GA>TG","GA>CG","GA>AG","GA>TC","GA>CC","GA>AC") # GA
    mutation_categories_complements <- c(mutation_categories_complements,"CA>TT","CA>GT","CA>AT","CA>TG","CA>GG","CA>AG","CA>TC","CA>GC","CA>AC") # CA
    mutation_categories_complements <- c(mutation_categories_complements,"AA>TT","AA>GT","AA>CT","AA>TG","AA>GG","AA>CG","AA>TC","AA>GC","AA>CC") # AA

    # identify doublets motif
    data$cat <- paste0(data$ref,">",data$alt)
    for(i in 1:length(mutation_categories_complements)) {
        data$cat[which(data$cat==mutation_categories_complements[i])] <- mutation_categories[i]
    }
    
    # count number of mutations per sample for each category
    mutation_categories <- data.table(cat=mutation_categories)
    data <- merge(mutation_categories[,list(cat)],data.table(sample=data$sample,cat=data$cat)[,.N,by=list(sample,cat)],by="cat",all=TRUE)
    data <- dcast(data,sample~cat,value.var="N")
    data <- data[!is.na(sample),drop=FALSE]
    data[is.na(data)] <- 0

    # make multi nucleotide counts matrix
    samples_names <- data$sample
    data <- as.matrix(data[,2:ncol(data),drop=FALSE])
    rownames(data) <- samples_names
    data <- data[sort(rownames(data)),,drop=FALSE]
    data <- data[,sort(colnames(data)),drop=FALSE]
    multi_nucleotides_counts <- array(0,c(nrow(data),78))
    rownames(multi_nucleotides_counts) <- rownames(data)
    colnames(multi_nucleotides_counts) <- sort(mutation_categories$cat)
    rows_contexts <- rownames(data)
    cols_contexts <- colnames(multi_nucleotides_counts)[which(colnames(multi_nucleotides_counts)%in%colnames(data))]
    multi_nucleotides_counts[rows_contexts,cols_contexts] <- data[rows_contexts,cols_contexts]

    # return multi nucleotides counts matrix
    return(multi_nucleotides_counts)

}
