# 1 factorial ANOVA over factor1
onefacaov_fac1 <- function(normRootdetectionTable, file_base = "1fac_ANOVA_factor1", draw_out = T) {
    
    fac2 <- levels(as.factor(normRootdetectionTable$Factor2))
    fac1 <- levels(as.factor(normRootdetectionTable$Factor1))
    nr_fac2 <- length(fac2)
    nr_fac1 <- length(fac1)
    matlist <- as.list(rep(NA, nr_fac2))
    
    for (tp in 1:nr_fac2) {
        # create subdataset with only one level of Factor2
        data.aov = normRootdetectionTable[which(normRootdetectionTable$Factor2 == fac2[tp]), ]
        # ANOVA
        lm.aov = lm(LengthMM ~ Factor1, data.aov, na.action = na.omit)
        pval.aov = anova(aov(lm.aov))[1, 5]
        # Tukey
        tuk = as.data.frame(TukeyHSD(aov(lm.aov), ordered = FALSE)$Factor1)
        tuk = cbind(rownames(tuk), tuk)
        colnames(tuk)[1] = "comparison"
        # get tukey comparison and put them in matrix
        mat = matrix(NA, nrow = nr_fac1, ncol = nr_fac1)
        colnames(mat) = fac1
        rownames(mat) = fac1
        rn = rownames(tuk)
        for (j in 1:(nr_fac1 - 1)) {
            for (k in (j + 1):nr_fac1) {
                # 
                idx = which(paste(fac1[k], "-", fac1[j], sep = "") == rn)
                # 
                if (length(idx) != 0) {
                  mat[j, k] = tuk[idx, 5]
                }
            }
        }
        matlist[[tp]] <- mat
        
        if (draw_out) {
            pdf(file = paste(file_base, "_", fac2[tp], ".pdf", sep = ""))
            # Draw the Tuky-p-value-Matrix
            col = matrix("black", nrow = nr_fac1, ncol = nr_fac1)
            col[lower.tri(col, diag = TRUE)] = "white"
            image(t(mat[nr_fac1:1, ]), col = c("red", "white"), breaks = c(0, 0.05, 1), axes = FALSE)
            title(main = paste(fac2[tp], ", pval.aov = ", round(pval.aov, 5), sep = ""))
            s = seq(from = 0, to = 1, length = nr_fac1)
            axis(1, at = s, labels = fac1, las = 2)
            axis(2, at = s, labels = rev(fac1), las = 2)
            abline(h = c(s - (s[2] - s[1])/2, s[nr_fac1] + (s[2] - s[1])/2), v = c(s - (s[2] - s[1])/2, s[nr_fac1] + (s[2] - s[1])/2), col = "grey")
            sg = expand.grid(rev(s), s)
            text(sg[, 2], sg[, 1], format(c(as.vector(mat), 0.123456789), digits = 2, nsmall = 3, scientiffic = FALSE)[1:nr_fac1^2], col = as.vector(col), cex = 0.5)
            dev.off()
        }
    }
    
    names(matlist) <- unique(normRootdetectionTable$Factor2)
    return(matlist)
}

# 1 factorial ANOVA over factor2
onefacaov_fac2 <- function(normRootdetectionTable, control, file_base = "1fac_ANOVA_factor2", draw_out = T) {
    
    # only last mat is returned --> list of dfs!!!
    
    fac2 <- levels(as.factor(normRootdetectionTable$Factor2))
    fac1 <- levels(as.factor(normRootdetectionTable$Factor1))
    nr_fac2 <- length(fac2)
    nr_fac1 <- length(fac1)
    con_position <- which(unique(normRootdetectionTable$Factor2) == control)
    treat_position <- which(unique(normRootdetectionTable$Factor2) != control)
    
    # create empty list to put in mats for every pair
    matl <- list()
    # loop over treat positions
    for (tp in treat_position) {
        mat = matrix(NA, nrow = 1, ncol = nr_fac1)
        rownames(mat) = c("pval")
        colnames(mat) = fac1
        
        # loop over factor1
        for (i in 1:nr_fac1) {
            # for each factor1: get control and one single treatment
            d = normRootdetectionTable[which(normRootdetectionTable$Factor1 == fac1[i] & (normRootdetectionTable$Factor2 == fac2[con_position] | normRootdetectionTable$Factor2 == 
                fac2[tp])), ]
            # ANOVA
            lm.all = lm(LengthMM ~ Factor2, d, na.action = na.omit)
            aov.all = aov(lm.all)
            mat[1, i] = anova(aov.all)[1, 5]
            
            # add mat to list with name
            name <- fac2[tp]
            matl[[name]] <- mat
            
            if (draw_out) {
                # Visualize the p-values
                col = matrix("black", nrow = 1, ncol = nr_fac1)
                col[lower.tri(col, diag = TRUE)] = "white"
                # filename
                pdf(file = paste(file_base, "_", fac2[tp], ".pdf", sep = ""), width = 6, height = 2.5)
                image(t(mat), col = c("red", "white"), breaks = c(0, 0.05, 1), axes = FALSE)
                title(main = fac2[tp])
                s = seq(from = 0, to = 1, length = nr_fac1)
                axis(1, at = s, labels = fac1, las = 2)
                axis(2, at = s[1], labels = "pval", las = 2)
                abline(v = c(s - (s[2] - s[1])/2, s[nr_fac1] + (s[2] - s[1])/2), col = "grey")
                sg = expand.grid(rev(s), s)
                text(s, s[1], round(as.vector(mat), 3), cex = 0.5)
                dev.off()
            }
        }
    }
    return(matl)
}

# aov all vs all
twofacaov <- function(normRootdetectionTable, file = "2fac_ANOVA_all_vs_all.pdf", draw_out = T, label_delim = ";") {
    
    # model that compares all vs all
    aov.all_vs_all <- aov(LengthMM ~ Factor1 * Factor2, data = normRootdetectionTable)
    # tukey post hoc test
    tuk = as.data.frame(TukeyHSD(aov.all_vs_all, ordered = FALSE)$`Factor1:Factor2`)
    tuk = cbind(rownames(tuk), tuk)
    # which and how many labels are there
    labs <- levels(normRootdetectionTable$Label)
    nr_labs <- length(labs)
    # create empty matrix
    mat = matrix(NA, nrow = nr_labs, ncol = nr_labs)
    colnames(mat) = labs
    rownames(mat) = labs
    rn = rownames(tuk)
    for (j in 1:(nr_labs - 1)) {
        for (k in (j + 1):nr_labs) {
            # get p-values and put them into the matrix
            idx = which(paste(gsub(label_delim, ":", labs[k]), "-", gsub(label_delim, ":", labs[j]), sep = "") == rn)
            if (length(idx) == 0) {
                idx = which(paste(gsub(label_delim, ":", labs[j]), "-", gsub(label_delim, ":", labs[k]), sep = "") == rn)
            }
            # 
            if (length(idx) != 0) {
                mat[j, k] = tuk[idx, 5]
            }
        }
    }
    if (draw_out) {
        # Draw the Tuky-p-value-Matrix
        col = matrix("black", nrow = nr_labs, ncol = nr_labs)
        col[lower.tri(col, diag = TRUE)] = "white"
        pdf(file = file)
        image(t(mat[nr_labs:1, ]), col = c("red", "white"), breaks = c(0, 0.05, 1), axes = FALSE)
        title(main = "ANOVA all vs all")
        s = seq(from = 0, to = 1, length = nr_labs)
        axis(1, at = s, labels = labs, las = 2, cex.axis = 0.6)
        axis(2, at = s, labels = rev(labs), las = 2, cex.axis = 0.6)
        abline(h = c(s - (s[2] - s[1])/2, s[nr_labs] + (s[2] - s[1])/2), v = c(s - (s[2] - s[1])/2, s[nr_labs] + (s[2] - s[1])/2), col = "grey")
        sg = expand.grid(rev(s), s)
        text(sg[, 2], sg[, 1], format(c(as.vector(mat), 0.123456789), digits = 2, nsmall = 3, scientiffic = FALSE)[1:nr_labs^2], col = as.vector(col), cex = 0.8)
        dev.off()
    }
    return(mat)
}

pairwise.2facaov <- function(normRootdetectionTable, control = "20", draw_out = T, file_base = "2fac_ANOVA_BH_corrected_", label_delim = ";") {
    
    normRootdetectionTable$Factor1 <- as.factor(normRootdetectionTable$Factor1)
    normRootdetectionTable$Factor2 <- as.factor(normRootdetectionTable$Factor2)
    
    
    # get all levels of factor2
    fac2 <- levels(as.factor(normRootdetectionTable$Factor2))
    # get level position of control and treatment of Factor2
    nr_fac2 <- length(fac2)
    con_position <- which(unique(normRootdetectionTable$Factor2) == control)
    treat_position <- which(unique(normRootdetectionTable$Factor2) != control)
    # get all levels of Factor1
    fac1 <- levels(normRootdetectionTable$Factor1)
    # hoe many levels are there in Factor1
    nr_fac1 <- length(fac1)
    
    
    dfl <- list()
    
    for (tp in treat_position) {
        # create empty matrix filled with NA / as factor1 names to rows and colums
        mat = matrix(NA, nrow = nr_fac1, ncol = nr_fac1)
        rownames(mat) = fac1
        colnames(mat) = fac1
        # loop over all pairs of Factor1 i = 1 - length(Factor1-1)
        
        for (i in 1:(nr_fac1 - 1)) {
            # j = i+1 - length(Factor1)
            for (j in (i + 1):nr_fac1) {
                # create data.frame with 2 ecotypes | 1. control 2. from loops
                d = normRootdetectionTable[which((normRootdetectionTable$Factor1 == fac1[i] | normRootdetectionTable$Factor1 == fac1[j]) & (normRootdetectionTable$Factor2 == 
                  fac2[con_position] | normRootdetectionTable$Factor2 == fac2[tp])), ]
                # definiere linear model for ANOVA
                lm.all = lm(log2(LengthMM) ~ Factor1 * Factor2, d, na.action = na.omit)
                # get matrix with all the p-values!!
                mat[i, j] = anova(aov(lm.all))[3, 5]
            }
        }
        
        # adjust p-values
        pvals = mat[upper.tri(mat)]
        mat.adj = matrix(NA, nrow = nr_fac1, ncol = nr_fac1)
        rownames(mat.adj) = fac1
        colnames(mat.adj) = fac1
        adjp = mt.rawp2adjp(pvals, proc = "BH")
        idx = order(adjp$index)
        mat.adj[upper.tri(mat.adj)] = adjp$adjp[idx, "BH"]
        mat = mat.adj
        if (draw_out) {
            # Visualisiation of P-Value-Matrix font color
            col = matrix("black", nrow = nr_fac1, ncol = nr_fac1)
            col[lower.tri(col, diag = TRUE)] = "white"
            # filename for pdf
            pdf(file = paste(file_base, control, "_vs_", fac2[tp], ".pdf", sep = ""))
            # draw matrix draw white and red cells
            image(t(mat[nr_fac1:1, ]), col = c("red", "white"), breaks = c(0, 0.05, 1), axes = FALSE)
            title(main = paste("treatment effect ", control, " vs ", fac2[tp]))
            # axis label
            s = seq(from = 0, to = 1, length = nr_fac1)
            axis(1, at = s, labels = fac1, las = 1, cex.axis = 0.5)
            axis(2, at = s, labels = rev(fac1), las = 1, cex.axis = 0.5)
            # grid lines
            abline(h = c(s - (s[2] - s[1])/2, s[nr_fac1] + (s[2] - s[1])/2), v = c(s - (s[2] - s[1])/2, s[nr_fac1] + (s[2] - s[1])/2), col = "grey")
            # print p-Werte into it
            sg = expand.grid(rev(s), s)
            text(sg[, 2], sg[, 1], cex = 0.3, format(c(as.vector(mat), 0.123456789), digits = 1, nsmall = 3, scientiffic = FALSE)[1:nr_fac1^2], col = as.vector(col))
            dev.off()
        }
        
        # Generate dataframes with p-values for each condition (for each loop) and put it in list dfl
        
        # get names from row and col
        rowCol <- expand.grid(rownames(mat), colnames(mat))
        # take only the upper.tri
        labs <- rowCol[as.vector(upper.tri(mat)), ]
        # add condition to labels
        labs[, 1] <- paste(labs[, 1], unique(normRootdetectionTable$Factor2)[tp], sep = label_delim)
        labs[, 2] <- paste(labs[, 2], unique(normRootdetectionTable$Factor2)[tp], sep = label_delim)
        # bind with values
        df <- cbind(labs, mat[upper.tri(mat)])
        # order data frame
        df <- df[, c(2, 1, 3)]
        # combine cells to new cell called name
        df <- unite(df, names, c(Var2, Var1), sep = "-", remove = TRUE)
        colnames(df)[2] <- "p.value"
        
        name <- fac2[tp]
        dfl[[name]] <- df
    }
    return(dfl)
}
