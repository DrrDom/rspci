#' Read and preprocess data about fragments contributions from a text file
#'
#' @param file_name name of the input text file with fragments contributions.
#' @param sep separator between molecule and fragment names. Default is ###.
#' @param keep_models character vector with model names to keep. Only for those
#' models, if their number is more than one, consensus (average) contributions
#' will be calculated.
#' @return melted data.frame with added consensus (average) contributions
#' @details An input file of an old format had '#' as a separator.
#' @export
#' @examples
#' file_name <- system.file("extdata", "free-wilson_frag_contributions.txt", package = "rspci")
#' df <- load_data(file_name)
#' @importFrom dplyr %>%
load_data <- function(file_name, sep = "###", keep_models = NULL) {
  # init contribution names for further replacement
  contrib_names <- c("overall", "hydrophobic", "electrostatic", "hydrogen bonding", "dispersive")
  names(contrib_names) <- c("overall", "LOGP", "CHARGE", "HB", "REFRACTIVITY")
  # load data
  df <- as.data.frame(data.table::fread(file_name, sep = "\t", header = T, na.strings = "", autostart = 3))
  # keep only selected models
  if (!is.null(keep_models)) {
    df <- df[grepl(paste0("^", keep_models, "_", collapse = "|"), df[, 1]), ]
  }
  # prepare df
  rnames <- colnames(df)[-1]
  colnames(df) <- NULL
  rownames(df) <- df[,1]
  df <- as.data.frame(t(df[,-1]))
  # get mol and frag names
  n <- as.data.frame(do.call(rbind, strsplit(rnames, sep)), stringsAsFactors = FALSE)
  colnames(n) <- c("MolID", "FragID")
  # add unique fragment id
  n$FragUID <- 1:nrow(n)
  # count occurencies
  nm <- n %>%
    dplyr::group_by(FragID) %>%
    dplyr::summarise(M = length(unique(MolID)), N = length(MolID))
  df <- cbind(n, nm[match(n$FragID, nm$FragID), -1], df)
  # melt
  df <- reshape2::melt(df, id.vars = 1:5)
  df <- cbind(df[, 1:5],
              do.call(rbind, strsplit(as.character(df$variable), "_")),
              df$value)
  colnames(df)[6:8] <- c("Model", "Property", "Contribution")
  df[, c("Model", "Property")] <- sapply(df[, c("Model", "Property")], as.character)
  # replace known descriptor property names
  v <- contrib_names[df[, "Property"]]
  v[is.na(v)] <- df[is.na(v), "Property"]
  df[, "Property"] <- v
  # add average consensus
  if (length(unique(df$Model)) > 1) {
    cons_df <- split(df, df$Model)[[1]]
    cons_df$Model <- "consensus"
    cons_df$Contribution <- rowMeans(sapply(split(df, df$Model), "[[", "Contribution"))
    df <- dplyr::bind_rows(df, cons_df)
  }
  return(df)
}



#' Read and preprocess data about fragments contributions from a text file created with spci-ext Python scripts
#'
#' @param file_name name of the input text file with fragments contributions.
#' @param sep separator between molecule and fragment names. Default is ###.
#' @return data.frame in long format
#' @details An input file doesn't have a header and consists of two columns with names and values
#' @export
#' @importFrom dplyr %>%
#' @examples
#' \dontrun{
#'   df <- load_data_ext(file_name)
#' }
load_data_ext <- function(file_name, sep = "###") {
  df <- as.data.frame(data.table::fread(file_name, sep = "\t", header = F, na.strings = "", autostart = 3))
  n <- as.data.frame(do.call(rbind, strsplit(df[, 1], sep)), stringsAsFactors = FALSE)
  colnames(n) <- c("MolID", "FragID")
  # count occurencies
  nm <- n %>%
    dplyr::group_by(FragID) %>%
    dplyr::summarise(M = length(unique(MolID)), N = length(MolID))
  df <- cbind(n, nm[match(n$FragID, nm$FragID), -1], Model = sub("\\..*$", "", basename(file_name)), Property = "overall", Contribution = df[, 2])
  return(df)
}



#' Add full fragments names as the first column to the existed data.frame
#'
#' @param df data.frame load with load_data function.
#' @param addM whether or not to add the count of molecules contatning each fragment
#' @param addN whether or not to add the count of each fragment across the whole data set
#' @return data.frame with the first column containing full fragments names
#' @details Convenience function to use instead of get_full_names. Input data.frame has
#' to contain three columns names \code{FragID}, \code{M} and \code{N}
#' @export
#' @examples
#' file_name <- system.file("extdata", "free-wilson_frag_contributions.txt", package = "rspci")
#' df <- load_data(file_name)
#' df <- add_full_names(df)
add_full_names <- function(df, addM = TRUE, addN = TRUE) {
  return(cbind(full_name = get_full_names(df, addM, addN), df))
}



#' Filter out fragments which are outside of specified ranges
#'
#' @param df input data.frame.
#' @param minM,maxM,minN,maxN values of minimum and maximum valus for M and N columns in the input data.frame
#' @return data.frame.
#' @export
#' @examples
#' file_name <- system.file("extdata", "free-wilson_frag_contributions.txt", package = "rspci")
#' df <- load_data(file_name)
#' df <- filter_by_frags_count(df)
filter_by_frags_count <- function(df, minM = 10, minN = 10, maxM = Inf, maxN = Inf) {
  dplyr::filter(df, N >= minN & N <= maxN & M >= minM & M <= maxM)
}



#' Keep only models which are specified
#'
#' @param df input data.frame.
#' @param keep_model_names character vector of models names
#' @return data.frame.
#' @details The input data.frame has to contain \code{Model} column.
#' @export
#' @examples
#' file_name <- system.file("extdata", "free-wilson_frag_contributions.txt", package = "rspci")
#' df <- load_data(file_name)
#' df <- filter_by_model_names(df, c("svm", "pls", "rf"))
filter_by_model_names <- function(df, keep_model_names) {
  df[df$Model %in% keep_model_names, ]
}



#' Keep only properties which are specified
#'
#' @param df input data.frame.
#' @param keep_prop_names character vector of models names
#' @param remove_prop_names character vector of models names
#' @return data.frame.
#' @details The input data.frame has to contain \code{Property} column.
#' @export
#' @examples
#' file_name <- system.file("extdata", "free-wilson_frag_contributions.txt", package = "rspci")
#' df <- load_data(file_name)
#' df1 <- filter_by_prop_names(df, "overall")
#' df2 <- filter_by_prop_names(df, remove_prop_names = c("overall", "dispersive"))
filter_by_prop_names <- function(df, keep_prop_names = NULL, remove_prop_names = NULL) {
  if (!is.null(keep_prop_names)) {
    df <- df[df$Property %in% keep_prop_names, ]
  }
  if (!is.null(remove_prop_names)) {
    df <- df[!(df$Property %in% remove_prop_names), ]
  }
  return(df)
}



#' Reorder fragments names represented as factor levels based on a specified metric
#'
#' @param df input data.frame.
#' @param by_model character name of a selected model name in column
#' @param by_prop character name of a selected property
#' @param frag_col_name character name of column containing desired fragments names to reorder
#' @param FUN metric functions applied to order different fragments
#' @return data.frame.
#' @details This function can be applied only to the melted input data.frame.
#' @export
#' @importFrom dplyr %>%
#' @examples
#' file_name <- system.file("extdata", "free-wilson_frag_contributions.txt", package = "rspci")
#' df <- load_data(file_name)
#' df <- reorder_data(df, "consensus", "overall", frag_col_name = "FragID")
reorder_data <- function(df, by_model, by_prop, FUN = median, frag_col_name = "full_name") {
  sorted_levels <- eval(substitute(
    df[df$Model == by_model & df$Property == by_prop, ] %>%
      dplyr::group_by(frag) %>%
      dplyr::summarise(res = FUN(Contribution)) %>%
      dplyr::arrange(res),
    list(frag = as.name(frag_col_name))))
  df[, frag_col_name] <- factor(df[, frag_col_name],
                                levels = get_col(sorted_levels, 1))
  return(df)
}



#' Add statistical significance to fragments contributions
#'
#' @param df input data.frame.
#' @param FUN function which returns p.value within the list of named values.
#' By default wilcox.test is used, t.test can be used also
#' @return data.frame with added two columns names pvalue and ptext:
#' with numerical p.value and with string representation of p.value
#' (*** < 0.001, ** < 0.01, * < 0.05)
#' @export
#' @importFrom stats wilcox.test
#' @importFrom dplyr %>%
#' @examples
#' file_name <- system.file("extdata", "free-wilson_frag_contributions.txt", package = "rspci")
#' df <- load_data(file_name)
#' df <- add_signif(df)
add_signif <- function(df, FUN = wilcox.test) {
  s <- df %>%
    dplyr::group_by(FragID, Model, Property) %>%
    dplyr::summarise(pvalue = round(FUN(Contribution)$p.value, 6)) %>%
    dplyr::mutate(ptext = 3 - findInterval(pvalue, c(0.001, 0.01, 0.05)))
  # for single 0 contribution p.value is NA, to avoid errors NAs are replaced with 0s
  s$ptext <- sapply(replace_na(s$ptext), function(i) paste0(rep("*", i), collapse = ""))
  dplyr::left_join(df, s)
}



#' Plot fragments contributions
#'
#' @param df input data.frame.
#' @param frag_name_col name of a column with fragments names to plot
#' @param contrib_col name of a column with fragments contribution values
#' @param plot_type maybe \code{"boxplot"} or \code{"barplot"}
#' @param FUN aggragate function used in barplot (should return named list of outputs with numeric p.value item)
#' @param show_sign_text column name of text represented significance (usually number of asterisks). If \code{NULL} no text will be plot.
#' @param show_sep_lines boolean to control of showing grey lines to separate different fragments on a plot
#' @param flip boolean to control the orientation of a plot
#' @param x_labels_angle rorate axis x text labels. It works only if \code{flip} was set to TRUE.
#' @return ggplot object which can be further modified before printing
#' @export
#' @importFrom stats median
#' @importFrom ggplot2 aes aes_string coord_flip element_blank element_line
#' element_rect element_text facet_wrap geom_boxplot geom_text geom_vline
#' ggplot position_dodge stat_summary theme xlab
#' @examples
#' file_name <- system.file("extdata", "free-wilson_frag_contributions.txt", package = "rspci")
#' df <- load_data(file_name)
#' df <- add_full_names(df)
#' df <- reorder_data(df, "consensus", "overall")
#' df <- add_signif(df)
#' plot_contrib(df)
#' plot_contrib(df, plot_type = "barplot", flip = FALSE, show_sign_text = "ptext")
plot_contrib <- function(df, frag_name_col = "full_name", contrib_col = "Contribution",
                         plot_type = "boxplot", FUN = median, show_sign_text = NULL,
                         show_sep_lines = TRUE, flip = TRUE, x_labels_angle = 60) {

  # load and detach ggplot2 package if it was not loaded by user
  # loaded <- "package:ggplot2" %in% search()
  # if (!loaded) {
  #   if(!require(ggplot2)) {
  #     return()
  #   }
  # }

  # tryCatch({

    g <- ggplot(df, aes_string(x = frag_name_col, y = contrib_col))

    # base plot and representation of several Properties as fill
    if (plot_type == "barplot") {
      if (length(unique(df$Property)) > 1) {
        g <- g + stat_summary(fun.y = FUN, geom = "bar", position = "dodge", aes(fill = Property))
      } else {
        g <- g + stat_summary(fun.y = FUN, geom = "bar", position = "dodge", fill = "darkgrey")
      }
    } else {
      if (plot_type == "boxplot") {
        if (length(unique(df$Property)) > 1) {
          g <- g + geom_boxplot(aes(fill = Property), outlier.size = 1, position = position_dodge(width = 0.9))
        } else {
          g <- g + geom_boxplot(outlier.size = 1)
        }
      }
    }

    # add significance text
    if (!is.null(show_sign_text)) {
      if ((!flip & length(unique(df$Property)) > 1) | (flip & length(unique(df$Property)) == 1)) {
        angle <- 90
      } else {
        angle <- 0
      }
      if (plot_type == "barplot") {
        g <- g + stat_summary(fun.y = FUN, geom = "text", position = position_dodge(width = 0.9), angle = angle,
                              aes_string(group = "Property", label = show_sign_text))
      }
      if (plot_type == "boxplot") {
        g <- g + geom_text(aes_string(group = "Property", label = show_sign_text, y = max(df[, contrib_col])),
                           position = position_dodge(width = 0.9), angle = angle)
      }
    }

    # add grid lines for separation of different fragments
    if (show_sep_lines) {
      g <- g + geom_vline(xintercept = 1.5:length(unique(df[, frag_name_col])), color = "grey")
    }

    # rename x axis
    g <- g + xlab("Fragment")

    # flip coordinates
    if (flip) {
      g <- g + coord_flip()
    }

    # process several models as grid
    if (length(unique(df$Model)) > 1) {
      if (flip) {
        g <- g + facet_wrap(~ Model, nrow = 1)
      } else {
        g <- g + facet_wrap(~ Model, ncol = 1)
      }
    }

    # visualization theme
    g <- g +
      theme(axis.title.y = element_text(size = 15),
            panel.background = element_rect(colour = "white", fill = "white"),
            panel.grid.major.x = element_line(color = "darkgrey", linetype = 3),
            panel.grid.minor.x = element_blank(),
            axis.line = element_line(size = 0.8),
            legend.position = "top",
            legend.key = element_rect(fill = NA))

    # update visualization in the case of not flipping
    if (!flip) {
      g <- g + theme(axis.text.x = element_text(angle = x_labels_angle, hjust = 1, vjust = 1),
                     panel.grid.major.y = element_line(linetype = 3, color = "darkgrey"),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.major.x = element_blank())
    }

    return(g)

  # }, finally = {if (!loaded) detach(package:ggplot2)})

}



#' Clustering fragment contributions
#'
#' @param data vector of fragment contributions.
#' @param molids vector of molecule IDs corresponding to fragment contributions
#' @return Mclust model object or NULL if data has a single unique observation value (see details).
#' @details Mclust model is a gaussian mixture model based on integrated complete-
#' data likelihood optimization criterion. If all values  in data are equal
#'  or a single value provided, then NULL is returned.
#' @export
#' @importFrom mclust mclustICL Mclust mclustBIC
#' @examples
#' file_name <- system.file("extdata", "BBB_frag_contributions.txt", package = "rspci")
#' df <- load_data(file_name)
#' dx <- dplyr::filter(df, FragID == "OH (aliphatic)", Model == "consensus", Property == "overall")
#' m <- clust(dx$Contribution, dx$MolID)
clust <- function(data, molids = NULL) {
  if (!is.null(molids)) {
    names(data) <- molids
  }
  if (length(unique(data)) > 1) {
    icl <- mclustICL(data, modelNames = "V")
    return(Mclust(data, G = which.max(icl), modelNames = "V"))
  } else {return(NULL)}
}



#' Parameters of each gaussian (cluster)
#'
#' @param model mclust model object
#' @return data.frame with mean, variance and proportion for each gaussian (cluster)
#' @export
#' @examples
#' file_name <- system.file("extdata", "BBB_frag_contributions.txt", package = "rspci")
#' df <- load_data(file_name)
#' dx <- dplyr::filter(df, FragID == "OH (aliphatic)", Model == "consensus", Property == "overall")
#' m <- clust(dx$Contribution, dx$MolID)
#' par <- get_clust_params(m)
get_clust_params <- function(model) {
  data.frame("mean" = model$parameters$mean,
             "variance" = model$parameters$variance$sigmasq,
             "proportion"=model$parameters$pro,
             row.names = unique(model$classification))
}



#' Get IDs of molecules for each cluster
#'
#' @param model mclust model
#' @param uncertainty the maximum level of uncertainty for molecules
#' to belong to a cluster. Molecules with uncertainty higher than this
#' threshold will not appear in the resulting list
#' @return list of vectors containing IDs of molecules belonging to
#' each cluster if IDs were supplied to \link{clust}. Otherwise indices
#' will be returned.
#' @export
#' @examples
#' file_name <- system.file("extdata", "BBB_frag_contributions.txt", package = "rspci")
#' df <- load_data(file_name)
#' dx <- dplyr::filter(df, FragID == "OH (aliphatic)", Model == "consensus", Property == "overall")
#' m <- clust(dx$Contribution, dx$MolID)
#' get_mol_ids(m, uncert = 0.2)
get_mol_ids <- function(model, uncertainty = 1) {
  if (!is.null(rownames(model$data))) {
    ids <- rownames(model$data)[model$uncertainty <= uncertainty]
  } else {
    ids <- 1:nrow(model$data)
  }
  lapply(split(ids, model$classification[model$uncertainty <= uncertainty]), unique)
}



#' Build mclust models for multiple fragments

#' @param data input data.frame
#' @param fragnames column of the data.frame containing fragments names (i.e. FragID or full_name)
#' @param contrib_col_name name of a column with contribution values
#' @param mol_col_name name of a column with names (ids) of molecules
#' @return list containing mclust models for fragments contained in data.frame
#' @details If all contributions of a fragment are identical the model for that fragment
#' will not be built.
#' @export
#' @examples
#' file_name <- system.file("extdata", "BBB_frag_contributions.txt", package = "rspci")
#' df <- load_data(file_name)
#' df <- dplyr::filter(df, Model == "consensus", Property == "overall")
#' df <- add_full_names(df)
#' models <- clust_all(df, "full_name")
clust_all <- function(data, fragnames, contrib_col_name = "Contribution", mol_col_name = "MolID") {
  m <- lapply(split(data, data[[fragnames]]), function(df) {
    clust(df[[contrib_col_name]], df[[mol_col_name]])
  })
  return(m[!sapply(m, is.null)])
}



#' Get number of gaussians (clusters) in mclust model
#'
#' @param model mclust model object
#' @return number of gaussians (clusters)
#' @export
#' @examples
#' file_name <- system.file("extdata", "BBB_frag_contributions.txt", package = "rspci")
#' df <- load_data(file_name)
#' dx <- dplyr::filter(df, FragID == "OH (aliphatic)", Model == "consensus", Property == "overall")
#' m <- clust(dx$Contribution, dx$MolID)
#' num <- get_num_clust(m)
get_num_clust <- function(model) {
  length(unique(model$classification))
}



#' Plot mclust model
#'
#' @param model mclust model object
#' @param binwidth width of bins on histogram
#' @param title plot title
#' @param xlab label of x axis
#' @param ylab label of y axis
#' @return Plot with the histogram of occurence of fragment contributions,
#' its kernel density esitmate (dashed) and gaussians detected by a model
#' (color solid)
#' @export
#' @importFrom ggplot2 aes geom_histogram geom_density
#' theme labs element_rect element_line element_text stat_function
#' @examples
#' file_name <- system.file("extdata", "BBB_frag_contributions.txt", package = "rspci")
#' df <- load_data(file_name)
#' dx <- dplyr::filter(df, FragID == "OH (aliphatic)", Model == "consensus", Property == "overall")
#' m <- clust(dx$Contribution, dx$MolID)
#' plot_mclust(m)

plot_mclust <- function(model, binwidth = 0.1, title = NULL, xlab = "Contribution", ylab = "density"){
  k <- model$G
  col2<- rep(c("#CC79A7", "#D55E00", "#0072B2", "#F0E442", "#009E73", "#56B4E9", "#E69F00"), len = k)

  # lt <- 1:k
  plot_component <- function(x, mu, sigsq, pr) {
    pr * dnorm(x, mean = mu, sd = sigsq)
  }
  p <-  ggplot(as.data.frame(model$data), aes(V1)) +
    geom_histogram(binwidth = binwidth, aes(y = ..density..), fill = "white", col = "black") +  geom_density(lty = 2, lwd = 1) +
    labs(title = title, x = xlab, y = ylab) +
    theme(panel.background = element_rect(fill = "white"),
          axis.line = element_line(color = "black"),
          plot.title = element_text(hjust = 0.5, size = 11))

  for (i in 1:k) {
    lt <- 1:k
    p <- p +
      stat_function(geom = "line",
                    fun =  plot_component,
                    args = list(mu = model$parameters$mean[i],
                                sigsq = sqrt(model$parameters$variance$sigmasq[i]),
                                pr = model$parameters$pro[i]), col = col2[i],
                    lwd = 1.5)
  }

  return (p)

}



#' Save plots of multiple mclust models in a grid image
#'
#' @param filename file name to save plots
#' @param models list of Mclust models
#' @param xlab label for x axis
#' @param ylab label for y axis
#' @export
#' @details If list contains more than 60 models, the function will return more
#'  than one png file, each next file will contain next (up to) 60 images. (File names will
#'  get suffices "_i", i = 1,2,3,...)
#' @importFrom ggplot2 ggsave
#' @importFrom gridExtra marrangeGrob
#' @examples
#' file_name <- system.file("extdata", "BBB_frag_contributions.txt", package = "rspci")
#' df <- load_data(file_name)
#' df <- dplyr::filter(df, Model == "consensus", Property == "overall")
#' df <- add_full_names(df)
#' models <- clust_all(df, "full_name")
#' save_mclust_plots("models.png", models)

save_mclust_plots <- function(filename, models, xlab = "contribution", ylab = "density") {

  pat <- rep(1:ceiling(length(models)/60), each = 60, length.out = length(models))

  if  (!is.null(names(models))) {
    fr_lst <- split(names(models), pat)
  } else {
    fr_lst <- split(seq(1:length(models)), pat)
  }
  if (length(fr_lst) == 1){
    plots <- lapply(fr_lst[[1]], function(nm) plot_mclust(models[[nm]], title = nm, xlab = NULL, ylab = NULL))
    ggsave(filename,
            marrangeGrob(grobs = plots,  ncol = 4, nrow = ceiling(length(plots)/4),
                         top = NULL, left=ylab, bottom=xlab),
            height  = 3.33* ceiling(length(plots)/4), width = 16)
  } else {
    count <- 1
    for (i in fr_lst) {

      plots <- lapply(i, function(nm) plot_mclust(models[[nm]], title = nm, xlab = NULL, ylab = NULL))

      ggsave(sub("^(.*)\\.png$", paste0("\\1_", count, ".png"), filename),
             marrangeGrob(grobs = plots,  ncol = 4, nrow = ceiling(length(plots)/4),
                          top = NULL, left=ylab, bottom=xlab),
             height  = 3.33* ceiling(length(plots)/4), width = 16)
      count <- count + 1

    }
  }

}


