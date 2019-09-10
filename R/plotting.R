
plot_hist <- function(normRootdetectionTable, draw_out = T, file = "data_distribution.pdf") {
    # other tests than shapiro-wilk
    
    plot_list <- list()
    
    # make sure Factor1 and Factor 2 are factors
    my_dataset_new$Factor1 <- as.factor(my_dataset_new$Factor1)
    my_dataset_new$Factor2 <- as.factor(my_dataset_new$Factor2)
    
    for (lev in 1:length(levels(my_dataset_new$Label))) {
        
        # create subset contain only the data for one label
        hist_sub <- subset(my_dataset_new, Label == levels(my_dataset_new$Label)[lev])
        hist_sub$Factor1 <- droplevels(hist_sub$Factor1)
        # shapiro-wilk test (> 0.05 significant normally distributed = green else (< 0.05) = failure --> red)
        shap <- shapiro.test(hist_sub$LengthMM)
        if (shap$p.value >= 0.05) {
            
            brx <- pretty(range(hist_sub$LengthMM), n = nclass.Sturges(hist_sub$LengthMM), min.n = 1)
            
            p <- ggplot(hist_sub, aes_string(x = hist_sub$LengthMM)) + geom_histogram(color = "black", fill = "green", breaks = brx) + theme_classic() + scale_x_continuous("LengthMM") + 
                annotate("text", x = Inf, y = Inf, label = round(shap$p.value, digits = 4), hjust = 1, vjust = 1, size = 3) + ggtitle(paste(unique(hist_sub$Factor1), 
                unique(hist_sub$Factor2, sep = " ")))
            
        } else {
            
            brx <- pretty(range(hist_sub$LengthMM), n = nclass.Sturges(hist_sub$LengthMM), min.n = 1)
            
            p <- ggplot(hist_sub, aes_string(x = hist_sub$LengthMM)) + geom_histogram(color = "black", fill = "red", breaks = brx) + theme_classic() + scale_x_continuous("LengthMM") + 
                annotate("text", x = Inf, y = Inf, label = round(shap$p.value, digits = 4), hjust = 1, vjust = 1, size = 3) + ggtitle(paste(unique(hist_sub$Factor1), 
                unique(hist_sub$Factor2, sep = " ")))
        }
        plot_list[[lev]] <- p
    }
    
    if (draw_out == T) {
        pdf(file)
        grid.arrange(grobs = plot_list, ncol = 3, nrow = 4)
        dev.off()
    }
    return(plot_list)
}


plot_abs <- function(normRootdetectionTable, plot_significance = F, twofacov_output, label_delim = ";", type = "jitter", plot_colours, letter_height = 2) {
    
    
    if (plot_significance) {
        
        # get Letters for twofacaov output
        mat_names <- character()
        mat_values <- numeric()
        # loop over matrix and get names + values
        for (j in 1:(length(row.names(twofacov_output)) - 1)) {
            for (k in (j + 1):length(colnames(twofacov_output))) {
                v <- twofacov_output[j, k]
                t <- paste(row.names(twofacov_output)[j], colnames(twofacov_output)[k], sep = "-")
                mat_names <- c(mat_names, t)
                mat_values <- c(mat_values, v)
            }
        }
        
        # combine names + values
        names(mat_values) <- mat_names
        # get df with letters and replace : with label delim!!
        letters <- data.frame(multcompLetters(mat_values)["Letters"])
        letters$Label <- gsub(":", label_delim, rownames(letters))
        
        
        # letter 1/5 of highest value
        y_coord <- ddply(my_dataset_new, .(Label), summarize, y = fivenum(LengthMM)[5])
        y_coord$y <- y_coord$y + letter_height
        plot_letters <- merge(letters, y_coord, by = "Label")
        
        if (type == "jitter") {
            
            
            # Input von oben : Farben, Größen
            abs_plot <- ggplot() + geom_boxplot(data = normRootdetectionTable, aes(x = normRootdetectionTable$Label, y = normRootdetectionTable$LengthMM), outlier.shape = NA) + 
                geom_jitter(data = normRootdetectionTable, aes(x = normRootdetectionTable$Label, y = normRootdetectionTable$LengthMM, colour = as.factor(normRootdetectionTable$Factor2)), 
                  position = position_jitter(0.2, seed = 1)) + geom_text(data = plot_letters, aes(x = plot_letters$Label, y = plot_letters$y, label = plot_letters$Letters)) + 
                labs(colour = "condition") + theme_classic() + scale_y_continuous(name = "hypocotyl length [mm]") + scale_x_discrete(name = "") + theme(axis.text.x = element_text(angle = 45, 
                hjust = 1, vjust = 1), legend.title = element_text(size = 16), legend.text = element_text(size = 12)) + stat_n_text(data = normRootdetectionTable, aes(x = normRootdetectionTable$Label, 
                y = normRootdetectionTable$LengthMM), angle = 90, size = 3)
            
        } else if (type == "box") {
            
            abs_plot <- ggplot() + geom_boxplot(data = normRootdetectionTable, aes(x = normRootdetectionTable$Label, y = normRootdetectionTable$LengthMM, fill = as.factor(normRootdetectionTable$Factor2))) + 
                geom_text(data = plot_letters, aes(x = plot_letters$Label, y = plot_letters$y, label = plot_letters$Letters)) + labs(fill = "condition") + theme_classic() + 
                scale_y_continuous(name = "hypocotyl length [mm]") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + stat_n_text(data = normRootdetectionTable, 
                aes(x = normRootdetectionTable$Label, y = normRootdetectionTable$LengthMM), angle = 90, size = 3)
            
        } else {
            stop("For type only \"box\" or \"jitter\" are yet implemented")
        }
        
    } else {
        
        if (type == "jitter") {
            
            # Input von oben : Farben, Größen
            abs_plot <- ggplot() + geom_boxplot(data = normRootdetectionTable, aes(x = normRootdetectionTable$Label, y = normRootdetectionTable$LengthMM), outlier.shape = NA) + 
                geom_jitter(data = normRootdetectionTable, aes(x = normRootdetectionTable$Label, y = normRootdetectionTable$LengthMM, colour = as.factor(normRootdetectionTable$Factor2)), 
                  position = position_jitter(0.2, seed = 1)) + labs(colour = "condition") + theme_classic() + scale_y_continuous(name = "hypocotyl length [mm]") + scale_x_discrete(name = "") + 
                theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.title = element_text(size = 16), legend.text = element_text(size = 12)) + 
                stat_n_text(data = normRootdetectionTable, aes(x = normRootdetectionTable$Label, y = normRootdetectionTable$LengthMM), angle = 90, size = 3)
            
        } else if (type == "box") {
            
            abs_plot <- ggplot() + geom_boxplot(data = normRootdetectionTable, aes(x = normRootdetectionTable$Label, y = normRootdetectionTable$LengthMM, fill = as.factor(normRootdetectionTable$Factor2))) + 
                labs(fill = "condition") + theme_classic() + scale_y_continuous(name = "hypocotyl length [mm]") + theme(axis.text.x = element_text(angle = 45, hjust = 1, 
                vjust = 1)) + stat_n_text(data = normRootdetectionTable, aes(x = normRootdetectionTable$Label, y = normRootdetectionTable$LengthMM), angle = 90, size = 3)
            
        } else {
            stop("For type only \"box\" or \"jitter\" are yet implemented")
        }
        
    }
    
    if (!missing(plot_colours)) {
        
        abs_plot <- abs_plot + scale_fill_manual(values = plot_colours) + scale_colour_manual(values = plot_colours)
    }
    
    print(abs_plot)
}


plot_rel <- function(normRootdetectionTable, pairwise.2facaov.Output, control = "20", plot_significance = T, type = "jitter") {
    
    # first step --> create relative data inside loop type zuweisung funktioniert nicht --> ungenutzes Arguement ??? Häää???
    relRootdetectionTable <- rel_data(normRootdetectionTable, control = control)
    
    if (plot_significance == F) {
        
        if (type == "jitter") {
            
            relative_plot <- ggplot() + geom_boxplot(data = relRootdetectionTable, aes(x = relRootdetectionTable$Label, y = relRootdetectionTable$relative_value), outlier.shape = NA) + 
                geom_jitter(data = relRootdetectionTable, aes(x = relRootdetectionTable$Label, y = relRootdetectionTable$relative_value, colour = as.factor(relRootdetectionTable$Label)), 
                  position = position_jitter(0.2, seed = 1)) + labs(colour = "ecotypes") + theme_classic() + scale_y_continuous(name = "% growth") + scale_x_discrete(name = "") + 
                theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.title = element_text(size = 16), legend.text = element_text(size = 12))
            
        } else if (type == "box") {
            
            relative_plot <- ggplot() + geom_boxplot(data = relRootdetectionTable, aes(x = relRootdetectionTable$Label, y = relRootdetectionTable$relative_value, fill = as.factor(relRootdetectionTable$Label))) + 
                labs(fill = "ecotypes") + theme_classic() + scale_y_continuous(name = "% growth") + scale_x_discrete(name = "") + theme(axis.text.x = element_text(angle = 45, 
                hjust = 1, vjust = 1), legend.title = element_text(size = 16), legend.text = element_text(size = 12))
            
        } else {
            stop("For type only box and jitter are allowed")
        }
        
        
    } else if (plot_significance == T) {
        
        if (type == "jitter") {
            # loop over length elements of 2facaov output and create lables and plot for every element
            relative_plot <- list()
            for (i in 1:length(pairwise.2facaov.Output)) {
                
                # get labels for every table
                labs <- data.frame(multcompLetters(structure(pairwise.2facaov.Output[[i]][, 2], names = as.character(pairwise.2facaov.Output[[i]][, 1])))["Letters"])
                labs$Label <- rownames(labs)
                y_coord <- ddply(relRootdetectionTable, .(Label), summarize, y = fivenum(relative_value)[5] + 10)
                plot_letter <- merge(labs, y_coord, by = "Label")
                
                
                rel_sub <- subset(relRootdetectionTable, Factor2 == names(pairwise.2facaov.Output)[i])
                rel_sub$Factor2 <- as.factor(rel_sub$Factor2)
                rel_sub$Factor2 <- droplevels(rel_sub$Factor2)
                
                # use annotate instead of geom_text() to use within a loop
                relative_plot_temp <- ggplot() + geom_boxplot(data = rel_sub, aes_string(x = rel_sub$Label, y = rel_sub$relative_value), outlier.shape = NA) + geom_jitter(data = rel_sub, 
                  aes_string(x = rel_sub$Label, y = rel_sub$relative_value, colour = as.factor(rel_sub$Label)), position = position_jitter(0.2, seed = 1)) + labs(colour = "ecotypes") + 
                  theme_classic() + scale_y_continuous(name = "% growth") + scale_x_discrete(name = "") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
                  legend.title = element_text(size = 16), legend.text = element_text(size = 12)) + annotate("text", x = plot_letter$Label, y = plot_letter$y, label = plot_letter$Letters)
                
                relative_plot[[i]] <- relative_plot_temp
            }
            
        } else if (type == "box") {
            
            # loop over length elements of 2facaov output and create lables and plot for every element
            relative_plot <- list()
            for (i in 1:length(pairwise.2facaov.Output)) {
                
                # get labels for every table
                labs <- data.frame(multcompLetters(structure(pairwise.2facaov.Output[[i]][, 2], names = as.character(pairwise.2facaov.Output[[i]][, 1])))["Letters"])
                labs$Label <- rownames(labs)
                y_coord <- ddply(relRootdetectionTable, .(Label), summarize, y = fivenum(relative_value)[5] + 10)
                plot_letter <- merge(labs, y_coord, by = "Label")
                
                
                rel_sub <- subset(relRootdetectionTable, Factor2 == names(pairwise.2facaov.Output)[i])
                rel_sub$Factor2 <- as.factor(rel_sub$Factor2)
                rel_sub$Factor2 <- droplevels(rel_sub$Factor2)
                
                # use annotate instead of geom_text() to use within a loop
                relative_plot_temp <- ggplot() + geom_boxplot(data = rel_sub, aes_string(x = rel_sub$Label, y = rel_sub$relative_value, fill = rel_sub$Label)) + labs(fill = "ecotypes") + 
                  theme_classic() + scale_y_continuous(name = "% growth") + scale_x_discrete(name = "") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
                  legend.title = element_text(size = 16), legend.text = element_text(size = 12)) + annotate("text", x = plot_letter$Label, y = plot_letter$y, label = plot_letter$Letters)
                
                relative_plot[[i]] <- relative_plot_temp
            }
        } else {
            stop("For type only box and jitter are allowed")
        }
        
    } else {
        stop("plot_significance can only be TRUE or FALSE")
    }
    
    return(relative_plot)
}
