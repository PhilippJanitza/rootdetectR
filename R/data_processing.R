


is.RootdetectionOutput <- function(RootdetectionOutput) {
    
    # sammelkommentar um weitere Abfragen zu sammeln sind 0 Werte vorhanden?
    
    d.f_bool <- is.data.frame(RootdetectionOutput)
    if (!d.f_bool) {
        stop("RootdetectionOutput is not a data.frame", call. = FALSE)
    }
    cols <- all(c("Label", "LengthMM", "LengthPx") %in% colnames(RootdetectionOutput))
    if (!cols) {
        stop("data.frame is missing one or more of the important columns Label, LenghtPx and LengthMM - check colnames!", call. = FALSE)
    }
    pxnum <- is.numeric(RootdetectionOutput$LengthPx)
    if (!pxnum) {
        stop("The column LenghtPx needs to store numeric values", call. = FALSE)
    }
    # print messages als warnings??
    std <- "10mm" %in% RootdetectionOutput$Label
    if (!std) {
        warning("10mm standard is missing")
        mmnum <- is.numeric(RootdetectionOutput$LengthMM)
        pxmm <- all(RootdetectionOutput$LengthMM == RootdetectionOutput$LengthPx)
        if (mmnum == T && pxmm == F) {
            warning("Are the LengthMM values already normalized with a standard?")
        } else {
            stop("LengthMM seems to be not normalized with a standard!")
        }
    }
    
    hyph <- any(grepl("-", RootdetectionOutput$Label))
    if (hyph == T) {
        stop("Labels containing '-' are not allowed!", call. = TRUE)
    }
    
    labna <- any(is.na(RootdetectionOutput$Label))
    if (labna == T) {
        stop("Labels containing NA!", call. = TRUE)
    }
    
    pxna <- any(is.na(RootdetectionOutput$LengthPx))
    if (pxna == T) {
        stop("LengthPx containing NA!", call. = TRUE)
    }
    
    pxzero <- any(RootdetectionOutput$LengthPx == 0)
    if (pxzero == T) {
        warning("LengthPx containg one or more 0")
    }
    
    
    if (all(c(d.f_bool, cols, pxnum)) && !all(hyph, labna, pxna)) {
        return(TRUE)
    } else {
        stop("The present input object does not fulfill the RootdetectionOutput standard.", call. = FALSE)
    }
}

norm_10mm_standard <- function(RootdetectionOutput, sort = TRUE, label_delim = ";") {
    
    # features --> change name of 10mm to something else choose length that should be used as normalisation (other than 10mm)
    
    # calc 10mm
    standard_mm <- subset(RootdetectionOutput, Label == "10mm")
    standard_mm_mean <- mean(RootdetectionOutput$LengthPx)
    RootdetectionOutput$LengthMM <- RootdetectionOutput$LengthPx/(standard_mm_mean/10)
    # delete 10mm
    RootdetectionOutput <- subset(RootdetectionOutput, Label != "10mm")
    RootdetectionOutput$Label <- droplevels(RootdetectionOutput$Label)
    
    # if sort is TRUE order the columns and devide labels according to label_delim (standard = ';')
    if (sort == TRUE) {
        RootdetectionOutput <- separate(data = RootdetectionOutput, col = Label, into = c("Factor1", "Factor2"), sep = label_delim, remove = F)
        RootdetectionOutput <- RootdetectionOutput[, c("Nr", "Filename", "RootNr", "Label", "Factor1", "Factor2", "LengthPx", "LengthMM")]
        
    }
    
    # return data.frame
    return(RootdetectionOutput)
    
}

se <- function(x) sd(x, na.rm = T)/sqrt(length(na.omit(x)))

summary_stat <- function(normRootdetectionTable, print_n = T) {
    
    # wählen was man alles haben möchte??
    sum <- ddply(normRootdetectionTable, .(Factor1, Factor2), summarize, n = length(LengthMM), median = median(LengthMM), mean = mean(LengthMM), sd = sd(LengthMM), se = se(LengthMM))
    return(sum)
    
}

normality_test <- function(normRootdetectionTable, test = "sw") {
    
    norm_table <- data.frame()
    
    for (lev in 1:length(levels(normRootdetectionTable$Label))) {
        # create subset contain only the data for one label
        hist_sub <- subset(normRootdetectionTable, Label == levels(normRootdetectionTable$Label)[lev])
        
        norm_table$Label <- c(norm_table$Label, levels(normRootdetectionTable$Label)[lev])
        norm_table$p.value <- c(norm_table$p.value, shapiro.test(hist_sub$LengthMM)$p.value)
        
    }
    
    return(shap_table)
}


rel_data <- function(normRootdetectionTable, control = "20") {
    # create subset containing only mock (control) data (name of the control was assigned in the beginning)
    rel_table_mock <- subset(normRootdetectionTable, Factor2 == control)
    # calculate median for all levels of Factor 1 and save to LengthMM_median_control
    rel_table_mock_median <- ddply(rel_table_mock, .(Factor1), summarize, LengthMM_median_control = median(LengthMM))
    # merge tables
    rel_table_merge <- merge(normRootdetectionTable, rel_table_mock_median, by = "Factor1")
    # new subset without the mock data
    rel_table <- subset(rel_table_merge, Factor2 != control)
    # calculate relative values for the new 'treatment-only' table
    rel_table$relative_value <- (100 * rel_table$LengthMM)/rel_table$LengthMM_median_control
    relRootdetectionTable <- rel_table[, c("Label", "Factor1", "Factor2", "RootNr", "relative_value")]
    
    return(relRootdetectionTable)
}
