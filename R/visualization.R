#' Plot observed Single Base Substitutions (SBS) counts for different groups of patients.
#'
#' @examples
#' data(plot_data_examples)
#' counts <- plot_data_examples[['groups.SBS.plot']][['counts']]
#' groups <- plot_data_examples[['groups.SBS.plot']][['groups']]
#' groupsSBSPlot(counts=counts,groups=groups)
#'
#' @title groupsSBSPlot
#' @param counts Matrix with Single Base Substitutions (SBS) counts data.
#' @param groups List where names are groups labels and elements are patients labels corresponding to rownames in counts.
#' @param normalize Boolean value; shall I normalize observed counts?
#' @param xlabels Boolean value; shall I display x labels?
#' @return A ggplot2 object.
#' @export groupsSBSPlot
#' @import ggplot2
#' @import gridExtra
#'
groupsSBSPlot <- function(counts, groups, normalize = TRUE,
    xlabels = FALSE) {

    # make the ggplot2 object
    glist <- list()
    for (i in seq_len(length(groups))) {
        plot_data <- counts[groups[[i]], ]
        if (normalize) {
            plot_data <- plot_data/rowSums(plot_data)
        }
        x_label <- NULL
        x_value <- NULL
        for (a in seq_len(nrow(plot_data))) {
            for (b in seq_len(ncol(plot_data))) {
                x_label <- c(x_label, colnames(plot_data)[b])
                x_value <- c(x_value, plot_data[a,
                  b])
            }
        }
        plot_data <- data.frame(Context = paste0(substr(x_label,
            1, 1), ".", substr(x_label, 7,
            7)), alt = paste0(substr(x_label,
            3, 3), ">", substr(x_label, 5,
            5)), value = x_value)
        plt <- ggplot(plot_data) + geom_boxplot(aes_string(x = "Context",
            y = "value", fill = "alt")) + facet_wrap(~alt,
            nrow = 1, scales = "free_x") +
            theme(axis.text.x = element_text(angle = 90,
                hjust = 1), panel.background = element_blank(),
                axis.line = element_line(colour = "black")) +
            ggtitle(names(groups)[i]) + theme(legend.position = "none") +
            ylab("Frequency of mutations")

        if (!xlabels) {
            plt <- plt + theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank())
        }

        glist[[i]] <- plt

    }

    # make the final plot
    grid.arrange(grobs = glist, ncol = ceiling(length(groups)/3))

}

#' Plot observed Multi-Nucleotide Variants (MNVs) counts for different groups of patients.
#'
#' @examples
#' data(plot_data_examples)
#' counts <- plot_data_examples[['groups.MNV.plot']][['counts']]
#' groups <- plot_data_examples[['groups.MNV.plot']][['groups']]
#' groupsMNVPlot(counts=counts,groups=groups)
#'
#' @title groupsMNVPlot
#' @param counts Matrix with Multi-Nucleotide Variants (MNVs) counts data.
#' @param groups List where names are groups labels and elements are patients labels corresponding to rownames in counts.
#' @param normalize Boolean value; shall I normalize observed counts?
#' @param xlabels Boolean value; shall I display x labels?
#' @return A ggplot2 object.
#' @export groupsMNVPlot
#' @import ggplot2
#' @import gridExtra
#'
groupsMNVPlot <- function(counts, groups, normalize = TRUE,
    xlabels = FALSE) {

    # make the ggplot2 object
    glist <- list()
    for (i in seq_len(length(groups))) {
        plot_data <- counts[groups[[i]], ]
        if (normalize) {
            plot_data <- plot_data/rowSums(plot_data)
        }
        x_label <- NULL
        x_value <- NULL
        for (a in seq_len(nrow(plot_data))) {
            for (b in seq_len(ncol(plot_data))) {
                x_label <- c(x_label, colnames(plot_data)[b])
                x_value <- c(x_value, plot_data[a,
                  b])
            }
        }
        plot_data <- data.frame(Context = substr(as.character(x_label),
            4, 5), alt = paste0(substr(as.character(x_label),
            1, 3), "NN"), value = x_value)
        plt <- ggplot(plot_data) + geom_boxplot(aes_string(x = "Context",
            y = "value", fill = "alt")) + facet_wrap(~alt,
            nrow = 1, scales = "free_x") +
            theme(axis.text.x = element_text(angle = 90,
                hjust = 1), panel.background = element_blank(),
                axis.line = element_line(colour = "black")) +
            ggtitle(names(groups)[i]) + theme(legend.position = "none") +
            ylab("Frequency of mutations")

        if (!xlabels) {
            plt <- plt + theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank())
        }

        glist[[i]] <- plt

    }

    # make the final plot
    grid.arrange(grobs = glist, ncol = ceiling(length(groups)/3))

}

#' Plot Single Base Substitutions (SBS) counts for a set of given patients.
#'
#' @examples
#' data(plot_data_examples)
#' counts <- plot_data_examples[['patients.SBS.plot']][['counts']]
#' patientsSBSPlot(trinucleotides_counts=counts,samples=rownames(counts)[seq_len(2)])
#'
#' @title patientsSBSPlot
#' @param trinucleotides_counts Trinucleotides counts matrix.
#' @param samples Name of the samples. This should match a rownames in trinucleotides_counts.
#' @param freq Boolean value; shall I display rates instead of counts?
#' @param xlabels Boolean value; shall I display x labels?
#' @return A ggplot2 object.
#' @export patientsSBSPlot
#' @import ggplot2
#' @import gridExtra
#' @importFrom data.table as.data.table :=
#' @importFrom reshape2 melt
#'
patientsSBSPlot <- function(trinucleotides_counts, samples = rownames(trinucleotides_counts),
    freq = FALSE, xlabels = FALSE) {

    # make samples data
    trinucleotides_counts <- trinucleotides_counts[samples, , drop = FALSE]
    if (freq) {
        trinucleotides_counts <- trinucleotides_counts/rowSums(trinucleotides_counts)
    }

    # separate context and alteration
    x <- as.data.table(melt(as.matrix(trinucleotides_counts), varnames = c("patient",
        "cat")))
    x[, `:=`("Context", paste0(substr(cat, 1, 1), ".", substr(cat, 7, 7)))]
    x[, `:=`("alt", paste0(substr(cat, 3, 3), ">", substr(cat, 5, 5)))]

    # make the ggplot2 object
    glist <- list()
    for (i in seq_len(nrow(trinucleotides_counts))) {
        plt <- ggplot(x[x$patient == rownames(trinucleotides_counts)[i]]) + geom_bar(aes_string(x = "Context",
            y = "value", fill = "alt"), stat = "identity", position = "identity") +
            facet_wrap(~alt, nrow = 1, scales = "free_x") + theme(axis.text.x = element_text(angle = 90,
            hjust = 1), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
            ggtitle(rownames(trinucleotides_counts)[i]) + theme(legend.position = "none") +
            ylab("Number of mutations")

        if (freq) {
            plt <- plt + ylab("Frequency of mutations")
        }

        if (!xlabels) {
            plt <- plt + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
        }

        glist[[i]] <- plt

    }

    # make the final plot
    grid.arrange(grobs = glist, ncol = ceiling(nrow(trinucleotides_counts)/3))

}

#' Plot Multi-Nucleotide Variants (MNVs) counts for a set of given patients.
#'
#' @examples
#' data(plot_data_examples)
#' counts <- plot_data_examples[['patients.MNV.plot']][['counts']]
#' patientsMNVPlot(multi_nucleotides_counts=counts,samples=rownames(counts)[seq_len(2)])
#'
#' @title patientsMNVPlot
#' @param multi_nucleotides_counts Multi-Nucleotide counts matrix.
#' @param samples Name of the samples. This should match a rownames in multi_nucleotides_counts
#' @param freq Boolean value; shall I display rates instead of counts?
#' @param xlabels Boolean value; shall I display x labels?
#' @return A ggplot2 object.
#' @export patientsMNVPlot
#' @import ggplot2
#' @import gridExtra
#' @importFrom data.table as.data.table :=
#' @importFrom reshape2 melt
#'
patientsMNVPlot <- function(multi_nucleotides_counts, samples = rownames(multi_nucleotides_counts),
    freq = FALSE, xlabels = FALSE) {

    # make samples data
    multi_nucleotides_counts <- multi_nucleotides_counts[samples, , drop = FALSE]
    if (freq) {
        multi_nucleotides_counts <- multi_nucleotides_counts/rowSums(multi_nucleotides_counts)
    }

    # separate context and alteration
    x <- as.data.table(melt(as.matrix(multi_nucleotides_counts), varnames = c("patient",
        "cat")))
    x[, `:=`("Context", substr(as.character(x$cat), 4, 5))]
    x[, `:=`("alt", paste0(substr(as.character(x$cat), 1, 3), "NN"))]

    # make the ggplot2 object
    glist <- list()
    for (i in seq_len(nrow(multi_nucleotides_counts))) {
        plt <- ggplot(x[x$patient == rownames(multi_nucleotides_counts)[i]]) + geom_bar(aes_string(x = "Context",
            y = "value", fill = "alt"), stat = "identity", position = "identity") +
            facet_wrap(~alt, nrow = 1, scales = "free_x") + theme(axis.text.x = element_text(angle = 90,
            hjust = 1), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
            ggtitle(rownames(multi_nucleotides_counts)[i]) + theme(legend.position = "none") +
            ylab("Number of mutations")

        if (freq) {
            plt <- plt + ylab("Frequency of mutations")
        }

        if (!xlabels) {
            plt <- plt + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
        }

        glist[[i]] <- plt

    }

    # make the final plot
    grid.arrange(grobs = glist, ncol = ceiling(nrow(multi_nucleotides_counts)/3))

}

#' Plot the inferred Single Base Substitutions (SBS) mutational signatures.
#'
#' @examples
#' data(plot_data_examples)
#' beta <- plot_data_examples[['signatures.SBS.plot']][['beta']]
#' signaturesSBSPlot(beta=beta)
#'
#' @title signaturesSBSPlot
#' @param beta Matrix with the inferred mutational signatures.
#' @param useRowNames Boolean value; shall I use the rownames from beta as names for the signatures?
#' @param xlabels Boolean value; shall I display x labels?
#' @return A ggplot2 object.
#' @export signaturesSBSPlot
#' @import ggplot2
#' @import gridExtra
#' @importFrom data.table as.data.table :=
#' @importFrom reshape2 melt
#'
signaturesSBSPlot <- function(beta, useRowNames = FALSE,
    xlabels = FALSE) {

    # set names of the signatures
    if (!useRowNames) {
        rownames(beta) <- paste0("Signature ", seq_len(nrow(beta)))
    }

    # separate context and alteration
    x <- as.data.table(melt(as.matrix(beta), varnames = c("signature",
        "cat")))
    x[, `:=`("Context", paste0(substr(cat, 1, 1), ".", substr(cat,
        7, 7)))]
    x[, `:=`("alt", paste0(substr(cat, 3, 3), ">", substr(cat,
        5, 5)))]

    # make the ggplot2 object
    glist <- list()
    for (i in seq_len(nrow(beta))) {
        plt <- ggplot(x[x$signature == rownames(beta)[i]]) +
            geom_bar(aes_string(x = "Context", y = "value",
                fill = "alt"), stat = "identity", position = "identity") +
            facet_wrap(~alt, nrow = 1, scales = "free_x") +
            theme(axis.text.x = element_text(angle = 90,
                hjust = 1), panel.background = element_blank(),
                axis.line = element_line(colour = "black")) +
            ggtitle(rownames(beta)[i]) + theme(legend.position = "none") +
            ylab("Frequency of mutations")

        if (!xlabels) {
            plt <- plt + theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank())
        }

        glist[[i]] <- plt

    }

    # make the final plot
    grid.arrange(grobs = glist, ncol = ceiling(nrow(beta)/3))

}

#' Plot the inferred Multi-Nucleotide Variants (MNVs) mutational signatures.
#'
#' @examples
#' data(plot_data_examples)
#' beta <- plot_data_examples[['signatures.MNV.plot']][['beta']]
#' signaturesMNVPlot(beta=beta)
#'
#' @title signaturesMNVPlot
#' @param beta Matrix with the inferred mutational signatures.
#' @param useRowNames Boolean value; shall I use the rownames from beta as names for the signatures?
#' @param xlabels Boolean value; shall I display x labels?
#' @return A ggplot2 object.
#' @export signaturesMNVPlot
#' @import ggplot2
#' @import gridExtra
#' @importFrom data.table as.data.table :=
#' @importFrom reshape2 melt
#'
signaturesMNVPlot <- function(beta, useRowNames = FALSE,
    xlabels = FALSE) {

    # set names of the signatures
    if (!useRowNames) {
        rownames(beta) <- paste0("Signature ", seq_len(nrow(beta)))
    }

    # separate context and alteration
    x <- as.data.table(melt(as.matrix(beta), varnames = c("signature",
        "cat")))
    x[, `:=`("Context", substr(as.character(x$cat), 4, 5))]
    x[, `:=`("alt", paste0(substr(as.character(x$cat), 1,
        3), "NN"))]

    # make the ggplot2 object
    glist <- list()
    for (i in seq_len(nrow(beta))) {
        plt <- ggplot(x[x$signature == rownames(beta)[i]]) +
            geom_bar(aes_string(x = "Context", y = "value",
                fill = "alt"), stat = "identity", position = "identity") +
            facet_wrap(~alt, nrow = 1, scales = "free_x") +
            theme(axis.text.x = element_text(angle = 90,
                hjust = 1), panel.background = element_blank(),
                axis.line = element_line(colour = "black")) +
            ggtitle(rownames(beta)[i]) + theme(legend.position = "none") +
            ylab("Frequency of mutations")

        if (!xlabels) {
            plt <- plt + theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank())
        }

        glist[[i]] <- plt

    }

    # make the final plot
    grid.arrange(grobs = glist, ncol = ceiling(nrow(beta)/3))

}

#' Plot observed Copy Number (CN) counts for different groups of patients.
#'
#' @examples
#' data(plot_data_examples)
#' counts <- plot_data_examples[['groups.CN.plot']][['counts']]
#' groups <- plot_data_examples[['groups.CN.plot']][['groups']]
#' groupsCNPlot(counts=counts,groups=groups)
#'
#' @title groupsCNPlot
#' @param counts Matrix with Copy Number (CN) counts data.
#' @param groups List where names are groups labels and elements are patients labels corresponding to rownames in counts.
#' @param normalize Boolean value; shall I normalize observed counts?
#' @param xlabels Boolean value; shall I display x labels?
#' @return A ggplot2 object.
#' @export groupsCNPlot
#' @import ggplot2
#' @import gridExtra
#'
groupsCNPlot <- function(counts, groups, normalize = TRUE,
    xlabels = FALSE) {

    # make the ggplot2 object
    glist <- list()
    for (i in seq_len(length(groups))) {
        plot_data <- counts[groups[[i]], ]
        if (normalize) {
            plot_data <- plot_data/rowSums(plot_data)
        }
        x_label <- NULL
        x_value <- NULL
        for (a in seq_len(nrow(plot_data))) {
            for (b in seq_len(ncol(plot_data))) {
                x_label <- c(x_label, colnames(plot_data)[b])
                x_value <- c(x_value, plot_data[a,
                  b])
            }
        }
        Context <- NULL
        for (j in as.character(x_label)) {
            Context <- c(Context, strsplit(j,
                ":")[[1]][[3]])
        }
        Alt <- NULL
        for (j in as.character(x_label)) {
            Alt <- c(Alt, paste0(toupper(strsplit(j,
                ":")[[1]][[2]]), " ", strsplit(j,
                ":")[[1]][[1]]))
        }
        Alt <- gsub("HOMDEL", "HD", Alt)
        Alt <- gsub("HET", "Het", Alt)
        plot_data <- data.frame(Context = factor(Context,
            levels = c("0-100kb", "100kb-1Mb",
                "1Mb-10Mb", "10Mb-40Mb", ">1Mb",
                ">40Mb")), alt = factor(Alt,
            levels = c("HD 0", "LOH 1", "LOH 2",
                "LOH 3-4", "LOH 5-8", "LOH 9+",
                "Het 2", "Het 3-4", "Het 5-8",
                "Het 9+")), value = x_value)
        plt <- ggplot(plot_data) + geom_boxplot(aes_string(x = "Context",
            y = "value", fill = "alt")) + facet_wrap(~alt,
            nrow = 1, scales = "free_x") +
            theme(axis.text.x = element_text(angle = 90,
                hjust = 1), panel.background = element_blank(),
                axis.line = element_line(colour = "black")) +
            ggtitle(names(groups)[i]) + theme(legend.position = "none") +
            ylab("Frequency of mutations")

        if (!xlabels) {
            plt <- plt + theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank())
        }

        glist[[i]] <- plt

    }

    # make the final plot
    grid.arrange(grobs = glist, ncol = ceiling(length(groups)/3))

}

#' Plot Copy Number (CN) counts for a set of given patients.
#'
#' @examples
#' data(plot_data_examples)
#' counts <- plot_data_examples[['patients.CN.plot']][['counts']]
#' patientsCNPlot(cn_data_counts=counts,samples=rownames(counts)[seq_len(2)])
#'
#' @title patientsCNPlot
#' @param cn_data_counts Copy Number counts matrix.
#' @param samples Name of the samples. This should match a rownames in cn_data_counts
#' @param freq Boolean value; shall I display rates instead of counts?
#' @param xlabels Boolean value; shall I display x labels?
#' @return A ggplot2 object.
#' @export patientsCNPlot
#' @import ggplot2
#' @import gridExtra
#' @importFrom data.table as.data.table :=
#' @importFrom reshape2 melt
#'
patientsCNPlot <- function(cn_data_counts, samples = rownames(cn_data_counts),
    freq = FALSE, xlabels = FALSE) {

    # make samples data
    cn_data_counts <- cn_data_counts[samples, , drop = FALSE]
    if (freq) {
        cn_data_counts <- cn_data_counts/rowSums(cn_data_counts)
    }

    # separate context and alteration
    x <- as.data.table(melt(as.matrix(cn_data_counts),
        varnames = c("patient", "cat")))
    Context <- NULL
    for (i in as.character(x$cat)) {
        Context <- c(Context, strsplit(i, ":")[[1]][[3]])
    }
    x[, `:=`("Context", factor(Context, levels = c("0-100kb",
        "100kb-1Mb", "1Mb-10Mb", "10Mb-40Mb", ">1Mb",
        ">40Mb")))]
    Alt <- NULL
    for (i in as.character(x$cat)) {
        Alt <- c(Alt, paste0(toupper(strsplit(i, ":")[[1]][[2]]),
            " ", strsplit(i, ":")[[1]][[1]]))
    }
    Alt <- gsub("HOMDEL", "HD", Alt)
    Alt <- gsub("HET", "Het", Alt)
    x[, `:=`("alt", factor(Alt, levels = c("HD 0",
        "LOH 1", "LOH 2", "LOH 3-4", "LOH 5-8", "LOH 9+",
        "Het 2", "Het 3-4", "Het 5-8", "Het 9+")))]

    # make the ggplot2 object
    glist <- list()
    for (i in seq_len(nrow(cn_data_counts))) {
        plt <- ggplot(x[x$patient == rownames(cn_data_counts)[i]]) +
            geom_bar(aes_string(x = "Context", y = "value",
                fill = "alt"), stat = "identity", position = "identity") +
            facet_wrap(~alt, nrow = 1, scales = "free_x") +
            theme(axis.text.x = element_text(angle = 90,
                hjust = 1), panel.background = element_blank(),
                axis.line = element_line(colour = "black")) +
            ggtitle(rownames(cn_data_counts)[i]) +
            theme(legend.position = "none") + ylab("Number of mutations")

        if (freq) {
            plt <- plt + ylab("Frequency of mutations")
        }

        if (!xlabels) {
            plt <- plt + theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank())
        }

        glist[[i]] <- plt

    }

    # make the final plot
    grid.arrange(grobs = glist, ncol = ceiling(nrow(cn_data_counts)/3))

}

#' Plot the inferred Copy Number (CN) mutational signatures.
#'
#' @examples
#' data(plot_data_examples)
#' beta <- plot_data_examples[['signatures.CN.plot']][['beta']]
#' signaturesCNPlot(beta=beta)
#'
#' @title signaturesCNPlot
#' @param beta Matrix with the inferred mutational signatures.
#' @param useRowNames Boolean value; shall I use the rownames from beta as names for the signatures?
#' @param xlabels Boolean value; shall I display x labels?
#' @return A ggplot2 object.
#' @export signaturesCNPlot
#' @import ggplot2
#' @import gridExtra
#' @importFrom data.table as.data.table :=
#' @importFrom reshape2 melt
#'
signaturesCNPlot <- function(beta, useRowNames = FALSE, xlabels = FALSE) {

    # set names of the signatures
    if (!useRowNames) {
        rownames(beta) <- paste0("Signature ", seq_len(nrow(beta)))
    }

    # separate context and alteration
    x <- as.data.table(melt(as.matrix(beta), varnames = c("signature",
        "cat")))
    Context <- NULL
    for (i in as.character(x$cat)) {
        Context <- c(Context, strsplit(i, ":")[[1]][[3]])
    }
    x[, `:=`("Context", factor(Context, levels = c("0-100kb",
        "100kb-1Mb", "1Mb-10Mb", "10Mb-40Mb", ">1Mb", ">40Mb")))]
    Alt <- NULL
    for (i in as.character(x$cat)) {
        Alt <- c(Alt, paste0(toupper(strsplit(i, ":")[[1]][[2]]),
            " ", strsplit(i, ":")[[1]][[1]]))
    }
    Alt <- gsub("HOMDEL", "HD", Alt)
    Alt <- gsub("HET", "Het", Alt)
    x[, `:=`("alt", factor(Alt, levels = c("HD 0", "LOH 1",
        "LOH 2", "LOH 3-4", "LOH 5-8", "LOH 9+", "Het 2",
        "Het 3-4", "Het 5-8", "Het 9+")))]

    # make the ggplot2 object
    glist <- list()
    for (i in seq_len(nrow(beta))) {
        plt <- ggplot(x[x$signature == rownames(beta)[i]]) +
            geom_bar(aes_string(x = "Context", y = "value",
                fill = "alt"), stat = "identity", position = "identity") +
            facet_wrap(~alt, nrow = 1, scales = "free_x") +
            theme(axis.text.x = element_text(angle = 90,
                hjust = 1), panel.background = element_blank(),
                axis.line = element_line(colour = "black")) +
            ggtitle(rownames(beta)[i]) + theme(legend.position = "none") +
            ylab("Frequency of mutations")

        if (!xlabels) {
            plt <- plt + theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank())
        }

        glist[[i]] <- plt

    }

    # make the final plot
    grid.arrange(grobs = glist, ncol = ceiling(nrow(beta)/3))

}

#' Plot observed Copy Number (Reduced, CX) counts for different groups of patients.
#'
#' @examples
#' data(plot_data_examples)
#' counts <- plot_data_examples[['groups.CX.plot']][['counts']]
#' groups <- plot_data_examples[['groups.CX.plot']][['groups']]
#' groupsCXPlot(counts=counts,groups=groups)
#'
#' @title groupsCXPlot
#' @param counts Matrix with Copy Number (Reduced, CX) counts data.
#' @param groups List where names are groups labels and elements are patients labels corresponding to rownames in counts.
#' @param normalize Boolean value; shall I normalize observed counts?
#' @param xlabels Boolean value; shall I display x labels?
#' @return A ggplot2 object.
#' @export groupsCXPlot
#' @import ggplot2
#' @import gridExtra
#'
groupsCXPlot <- function(counts, groups, normalize = TRUE,
    xlabels = FALSE) {

    # make the ggplot2 object
    Context_values <- c("S1", "S2", "S3", "S4",
        "S5", "S6", "S7", "S8", "S9", "S10",
        "S11", "S12", "S13", "S14", "S15",
        "S16", "S17", "S18", "S19", "S20",
        "S21", "S22", "C1", "C2", "C3", "C4",
        "C5", "C6", "C7", "C8", "C9", "C10",
        "10MB1", "10MB2", "10MB3", "CHRARM1",
        "CHRARM2", "CHRARM3", "CHRARM4", "CHRARM5",
        "CN1", "CN2", "CN3")
    Context <- factor(Context_values, levels = Context_values)
    Alt_values <- c("Segment size", "Segment size",
        "Segment size", "Segment size", "Segment size",
        "Segment size", "Segment size", "Segment size",
        "Segment size", "Segment size", "Segment size",
        "Segment size", "Segment size", "Segment size",
        "Segment size", "Segment size", "Segment size",
        "Segment size", "Segment size", "Segment size",
        "Segment size", "Segment size", "Changepoint",
        "Changepoint", "Changepoint", "Changepoint",
        "Changepoint", "Changepoint", "Changepoint",
        "Changepoint", "Changepoint", "Changepoint",
        "Breakpoints", "Breakpoints", "Breakpoints",
        "Breakpoints", "Breakpoints", "Breakpoints",
        "Breakpoints", "Breakpoints", "Oscillating",
        "Oscillating", "Oscillating")
    alt <- factor(Alt_values, levels = unique(Alt_values))
    colnames(counts) <- paste0(alt, " ", Context)
    glist <- list()
    for (i in seq_len(length(groups))) {
        plot_data <- counts[groups[[i]], ]
        if (normalize) {
            plot_data <- plot_data/rowSums(plot_data)
        }
        x_label <- NULL
        x_value <- NULL
        for (a in seq_len(nrow(plot_data))) {
            for (b in seq_len(ncol(plot_data))) {
                x_label <- c(x_label, colnames(plot_data)[b])
                x_value <- c(x_value, plot_data[a,
                  b])
            }
        }
        Context_values <- NULL
        Alt_values <- NULL
        for (j in as.character(x_label)) {
            if (strsplit(j, " ")[[1]][[1]] ==
                "Segment") {
                Context_values <- c(Context_values,
                  strsplit(j, " ")[[1]][[3]])
                Alt_values <- c(Alt_values,
                  "Segment size")
            } else {
                Context_values <- c(Context_values,
                  strsplit(j, " ")[[1]][[2]])
                Alt_values <- c(Alt_values,
                  strsplit(j, " ")[[1]][[1]])
            }
        }
        Context <- factor(Context_values, levels = c("S1",
            "S2", "S3", "S4", "S5", "S6", "S7",
            "S8", "S9", "S10", "S11", "S12",
            "S13", "S14", "S15", "S16", "S17",
            "S18", "S19", "S20", "S21", "S22",
            "C1", "C2", "C3", "C4", "C5", "C6",
            "C7", "C8", "C9", "C10", "10MB1",
            "10MB2", "10MB3", "CHRARM1", "CHRARM2",
            "CHRARM3", "CHRARM4", "CHRARM5",
            "CN1", "CN2", "CN3"))
        alt <- factor(Alt_values, levels = c("Segment size",
            "Changepoint", "Breakpoints", "Oscillating"))
        plot_data <- data.frame(Context = Context,
            alt = alt, value = x_value)
        plt <- ggplot(plot_data) + geom_boxplot(aes_string(x = "Context",
            y = "value", fill = "alt")) + facet_wrap(~alt,
            nrow = 1, scales = "free_x") +
            theme(axis.text.x = element_text(angle = 90,
                hjust = 1), panel.background = element_blank(),
                axis.line = element_line(colour = "black")) +
            ggtitle(names(groups)[i]) + theme(legend.position = "none") +
            ylab("Frequency of mutations")

        if (!xlabels) {
            plt <- plt + theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank())
        }

        glist[[i]] <- plt

    }

    # make the final plot
    grid.arrange(grobs = glist, ncol = ceiling(length(groups)/3))

}

#' Plot Copy Number (Reduced, CX) counts for a set of given patients.
#'
#' @examples
#' data(plot_data_examples)
#' counts <- plot_data_examples[['patients.CX.plot']][['counts']]
#' patientsCXPlot(cn_data_counts=counts,samples=rownames(counts)[seq_len(2)])
#'
#' @title patientsCXPlot
#' @param cn_data_counts Copy Number counts matrix.
#' @param samples Name of the samples. This should match a rownames in cn_data_counts
#' @param freq Boolean value; shall I display rates instead of counts?
#' @param xlabels Boolean value; shall I display x labels?
#' @return A ggplot2 object.
#' @export patientsCXPlot
#' @import ggplot2
#' @import gridExtra
#' @importFrom data.table as.data.table :=
#' @importFrom reshape2 melt
#'
patientsCXPlot <- function(cn_data_counts, samples = rownames(cn_data_counts),
    freq = FALSE, xlabels = FALSE) {

    # make samples data
    cn_data_counts <- cn_data_counts[samples, , drop = FALSE]
    if (freq) {
        cn_data_counts <- cn_data_counts/rowSums(cn_data_counts)
    }

    # separate context and alteration
    Context_values <- c("S1", "S2", "S3", "S4", "S5",
        "S6", "S7", "S8", "S9", "S10", "S11", "S12",
        "S13", "S14", "S15", "S16", "S17", "S18", "S19",
        "S20", "S21", "S22", "C1", "C2", "C3", "C4",
        "C5", "C6", "C7", "C8", "C9", "C10", "10MB1",
        "10MB2", "10MB3", "CHRARM1", "CHRARM2", "CHRARM3",
        "CHRARM4", "CHRARM5", "CN1", "CN2", "CN3")
    Context <- factor(Context_values, levels = Context_values)
    Alt_values <- c("Segment size", "Segment size",
        "Segment size", "Segment size", "Segment size",
        "Segment size", "Segment size", "Segment size",
        "Segment size", "Segment size", "Segment size",
        "Segment size", "Segment size", "Segment size",
        "Segment size", "Segment size", "Segment size",
        "Segment size", "Segment size", "Segment size",
        "Segment size", "Segment size", "Changepoint",
        "Changepoint", "Changepoint", "Changepoint",
        "Changepoint", "Changepoint", "Changepoint",
        "Changepoint", "Changepoint", "Changepoint",
        "Breakpoints", "Breakpoints", "Breakpoints",
        "Breakpoints", "Breakpoints", "Breakpoints",
        "Breakpoints", "Breakpoints", "Oscillating",
        "Oscillating", "Oscillating")
    alt <- factor(Alt_values, levels = unique(Alt_values))
    colnames(cn_data_counts) <- paste0(alt, " ", Context)
    x <- as.data.table(melt(as.matrix(cn_data_counts),
        varnames = c("patient", "cat")))
    Context_values <- NULL
    Alt_values <- NULL
    for (i in as.character(x$cat)) {
        if (strsplit(i, " ")[[1]][[1]] == "Segment") {
            Context_values <- c(Context_values, strsplit(i,
                " ")[[1]][[3]])
            Alt_values <- c(Alt_values, "Segment size")
        } else {
            Context_values <- c(Context_values, strsplit(i,
                " ")[[1]][[2]])
            Alt_values <- c(Alt_values, strsplit(i,
                " ")[[1]][[1]])
        }
    }
    Context <- factor(Context_values, levels = c("S1",
        "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9",
        "S10", "S11", "S12", "S13", "S14", "S15", "S16",
        "S17", "S18", "S19", "S20", "S21", "S22", "C1",
        "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9",
        "C10", "10MB1", "10MB2", "10MB3", "CHRARM1",
        "CHRARM2", "CHRARM3", "CHRARM4", "CHRARM5",
        "CN1", "CN2", "CN3"))
    alt <- factor(Alt_values, levels = c("Segment size",
        "Changepoint", "Breakpoints", "Oscillating"))
    x[, `:=`("Context", Context)]
    x[, `:=`("alt", alt)]

    # make the ggplot2 object
    glist <- list()
    for (i in seq_len(nrow(cn_data_counts))) {
        plt <- ggplot(x[x$patient == rownames(cn_data_counts)[i]]) +
            geom_bar(aes_string(x = "Context", y = "value",
                fill = "alt"), stat = "identity", position = "identity") +
            facet_wrap(~alt, nrow = 1, scales = "free_x") +
            theme(axis.text.x = element_text(angle = 90,
                hjust = 1), panel.background = element_blank(),
                axis.line = element_line(colour = "black")) +
            ggtitle(rownames(cn_data_counts)[i]) +
            theme(legend.position = "none") + ylab("Number of mutations")

        if (freq) {
            plt <- plt + ylab("Frequency of mutations")
        }

        if (!xlabels) {
            plt <- plt + theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank())
        }

        glist[[i]] <- plt

    }

    # make the final plot
    grid.arrange(grobs = glist, ncol = ceiling(nrow(cn_data_counts)/3))

}

#' Plot the inferred Copy Number (Reduced, CX) mutational signatures.
#'
#' @examples
#' data(plot_data_examples)
#' beta <- plot_data_examples[['signatures.CX.plot']][['beta']]
#' signaturesCXPlot(beta=beta)
#'
#' @title signaturesCXPlot
#' @param beta Matrix with the inferred mutational signatures.
#' @param useRowNames Boolean value; shall I use the rownames from beta as names for the signatures?
#' @param xlabels Boolean value; shall I display x labels?
#' @return A ggplot2 object.
#' @export signaturesCXPlot
#' @import ggplot2
#' @import gridExtra
#' @importFrom data.table as.data.table :=
#' @importFrom reshape2 melt
#'
signaturesCXPlot <- function(beta, useRowNames = FALSE, xlabels = FALSE) {

    # set names of the signatures
    if (!useRowNames) {
        rownames(beta) <- paste0("Signature ", seq_len(nrow(beta)))
    }

    # separate context and alteration
    Context_values <- c("S1", "S2", "S3", "S4", "S5", "S6",
        "S7", "S8", "S9", "S10", "S11", "S12", "S13", "S14",
        "S15", "S16", "S17", "S18", "S19", "S20", "S21",
        "S22", "C1", "C2", "C3", "C4", "C5", "C6", "C7",
        "C8", "C9", "C10", "10MB1", "10MB2", "10MB3", "CHRARM1",
        "CHRARM2", "CHRARM3", "CHRARM4", "CHRARM5", "CN1",
        "CN2", "CN3")
    Context <- factor(Context_values, levels = Context_values)
    Alt_values <- c("Segment size", "Segment size", "Segment size",
        "Segment size", "Segment size", "Segment size", "Segment size",
        "Segment size", "Segment size", "Segment size", "Segment size",
        "Segment size", "Segment size", "Segment size", "Segment size",
        "Segment size", "Segment size", "Segment size", "Segment size",
        "Segment size", "Segment size", "Segment size", "Changepoint",
        "Changepoint", "Changepoint", "Changepoint", "Changepoint",
        "Changepoint", "Changepoint", "Changepoint", "Changepoint",
        "Changepoint", "Breakpoints", "Breakpoints", "Breakpoints",
        "Breakpoints", "Breakpoints", "Breakpoints", "Breakpoints",
        "Breakpoints", "Oscillating", "Oscillating", "Oscillating")
    alt <- factor(Alt_values, levels = unique(Alt_values))
    colnames(beta) <- paste0(alt, " ", Context)
    x <- as.data.table(melt(as.matrix(beta), varnames = c("signature",
        "cat")))
    Context_values <- NULL
    Alt_values <- NULL
    for (i in as.character(x$cat)) {
        if (strsplit(i, " ")[[1]][[1]] == "Segment") {
            Context_values <- c(Context_values, strsplit(i,
                " ")[[1]][[3]])
            Alt_values <- c(Alt_values, "Segment size")
        } else {
            Context_values <- c(Context_values, strsplit(i,
                " ")[[1]][[2]])
            Alt_values <- c(Alt_values, strsplit(i, " ")[[1]][[1]])
        }
    }
    Context <- factor(Context_values, levels = c("S1", "S2",
        "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10",
        "S11", "S12", "S13", "S14", "S15", "S16", "S17",
        "S18", "S19", "S20", "S21", "S22", "C1", "C2", "C3",
        "C4", "C5", "C6", "C7", "C8", "C9", "C10", "10MB1",
        "10MB2", "10MB3", "CHRARM1", "CHRARM2", "CHRARM3",
        "CHRARM4", "CHRARM5", "CN1", "CN2", "CN3"))
    alt <- factor(Alt_values, levels = c("Segment size",
        "Changepoint", "Breakpoints", "Oscillating"))
    x[, `:=`("Context", Context)]
    x[, `:=`("alt", alt)]

    # make the ggplot2 object
    glist <- list()
    for (i in seq_len(nrow(beta))) {
        plt <- ggplot(x[x$signature == rownames(beta)[i]]) +
            geom_bar(aes_string(x = "Context", y = "value",
                fill = "alt"), stat = "identity", position = "identity") +
            facet_wrap(~alt, nrow = 1, scales = "free_x") +
            theme(axis.text.x = element_text(angle = 90,
                hjust = 1), panel.background = element_blank(),
                axis.line = element_line(colour = "black")) +
            ggtitle(rownames(beta)[i]) + theme(legend.position = "none") +
            ylab("Frequency of mutations")

        if (!xlabels) {
            plt <- plt + theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank())
        }

        glist[[i]] <- plt

    }

    # make the final plot
    grid.arrange(grobs = glist, ncol = ceiling(nrow(beta)/3))

}
