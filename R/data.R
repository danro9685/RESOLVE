#' @name patients
#' @title Point mutations for 560 breast tumors
#' @description Dataset of counts of the point mutations detected in 560 breast tumors published in Nik-Zainal, Serena, et al. (2016).
#' @docType data
#' @usage data(patients)
#' @format Counts of the point mutations
#' @source Nature. 2016 Jun 2;534(7605):47-54 (https://www.nature.com/articles/nature17676).
#' @return Counts of point mutations for 560 tumors and 96 trinucleotides
NULL

#' @name background
#' @title Germline replication error
#' @description Germline replication error estimated in Rahbari, Raheleh, et al. (2016).
#' @docType data
#' @usage data(background)
#' @format Vector of rates
#' @source Nat Genet. 2016 Feb;48(2):126-133 (https://www.nature.com/articles/ng.3469).
#' @return Vector of rates for the 96 trinucleotides
NULL

#' @name background2
#' @title COSMIC replication error
#' @description Background replication error signature derived from COSMIC SBS5.
#' @docType data
#' @usage data(background2)
#' @format Vector of rates
#' @source COSMIC database (https://cancer.sanger.ac.uk/cosmic/signatures) v3.
#' @return Vector of rates for the 96 trinucleotides
NULL

#' @name ssm560_reduced
#' @title A reduced version of the point mutations for 560 breast tumors in the format compatible with the import function
#' @description Reduced versione of the dataset of counts of the point mutations detected in 560 breast tumors published in Nik-Zainal, Serena, et al. (2016).
#' @docType data
#' @usage data(ssm560_reduced)
#' @format Reduced versione of the counts of the point mutations in the format compatible with the import function
#' @source Nature. 2016 Jun 2;534(7605):47-54 (https://www.nature.com/articles/nature17676).
#' @return Reduced versione of the counts of point mutations for 560 tumors and 96 trinucleotides in the format compatible with the import function
NULL

#' @name cn_example_reduced
#' @title A reduced version of the copy number data for 5 TCGA samples in the format compatible with the import function
#' @description Reduced versione of the dataset of counts of copy numbers detected in TCGA tumors published in Steele, Christopher D., et al. (2022).
#' @docType data
#' @usage data(cn_example_reduced)
#' @format Reduced versione of the counts of copy numbers in the format compatible with the import function
#' @source Nature. 2022 Jun;606(7916):984-991 (https://www.nature.com/articles/s41586-022-04738-6).
#' @return Reduced versione of the counts of copy numbers for 5 tumors and 48 copy number classes in the format compatible with the import function
NULL

#' @name plot_data_examples
#' @title List data structure to run examples
#' @description List data structure to run examples.
#' @docType data
#' @usage data(plot_data_examples)
#' @format List data structure to run examples
#' @source List data structure to run examples.
#' @return List data structure to run examples
NULL
