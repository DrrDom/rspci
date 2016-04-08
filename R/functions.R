#' Read and preprocess data about fragments contributions from a text file
#'
#' @param file_name name of the input text file with fragments contributions.
#' @param sep separator between molecule and fragment names. Default is ###.
#' @param keep_model character vector with model names to keep. Only for those
#' models, if their number is more than one, consensus (average) contributions
#' will be calculated.
#' @return melted data.frame with added consensus (average) contributions
#' @details An input file of an old format had '#' as a separator.
#' @export
#' @examples
#' file_name <- system.file("extdata", "free-wilson_frag_contributions.txt", package = "rspci")
#' df <- load_data(file_name)
load_data <- function(file_name, sep = "###", keep_models = NULL) {
  # init contribution names for further replacement
  contrib_names <- c("overall", "hydrophobic", "electrostatic", "hydrogen bonding", "dispersive")
  names(contrib_names) <- c("overall", "LOGP", "CHARGE", "HB", "REFRACTIVITY")
  # load data
  df <- as.data.frame(data.table::fread(file_name, sep = "\t", header = T, na.strings = "", autostart = 3))
  rnames <- colnames(df)[-1]
  colnames(df) <- NULL
  rownames(df) <- df[,1]
  df <- as.data.frame(t(df[,-1]))
  # get mol and frag names
  n <- as.data.frame(do.call(rbind, strsplit(rnames, sep)), stringsAsFactors = FALSE)
  colnames(n) <- c("MolID", "FragID")
  # count occurencies
  nm <- n %>%
    dplyr::group_by(FragID) %>%
    dplyr::summarise(M = length(unique(MolID)), N = length(MolID))
  df <- cbind(n, nm[match(n$FragID, nm$FragID), -1], df)
  # melt
  df <- reshape2::melt(df, id.vars = 1:4)
  df <- cbind(df[, 1:4],
              do.call(rbind, strsplit(as.character(df$variable), "_")),
              df$value)
  colnames(df)[5:7] <- c("Model", "Property", "Contribution")
  df[, c("Model", "Property")] <- sapply(df[, c("Model", "Property")], as.character)
  # replace known descriptor property names
  v <- contrib_names[df[, "Property"]]
  v[is.na(v)] <- df[is.na(v), "Property"]
  df[, "Property"] <- v
  # keep only selected models
  if (!is.null(keep_models)) {
    df <- df %>% dplyr::filter(Model %in% keep_models)
  }
  # add average consensus
  if (length(unique(df$Model)) > 1) {
    avg <- df %>%
      dplyr::group_by(MolID, FragID, M, N, Property) %>%
      dplyr::summarise(Contribution = mean(Contribution))
    df <- dplyr::bind_rows(df, data.frame(avg[, 1:4], Model = "consensus", avg[, 5:6]))
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
#' @examples
#' df <- load_data_ext(file_name)
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
#' @return data.frame with the first column containing full fragments names
#' @details Convenience function to use instead of get_full_names. Input data.frame has
#' to contain three columns names \code{FragID}, \code{M} and \code{N}
#' @export
#' @examples
#' file_name <- system.file("extdata", "free-wilson_frag_contributions.txt", package = "rspci")
#' df <- load_data(file_name)
#' df <- add_full_names(df)
add_full_names <- function(df, addM = TRUE, addN = TRUE) {
  return(cbind(full_name = rspci:::get_full_names(df, addM, addN), df))
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
#' @param by_value character name of a selected value
#' @param frag_name character name of column containing desired fragments names to reorder
#' @param FUN metric functions applied to order different fragments
#' @return data.frame.
#' @details This function can be applied only to the melted input data.frame.
#' @export
#' @examples
#' file_name <- system.file("extdata", "free-wilson_frag_contributions.txt", package = "rspci")
#' df <- load_data(file_name)
#' df <- reorder_data(df, frag_col_name = "FragID")
reorder_data <- function(df, by_model, by_prop, FUN = median, frag_col_name = "full_name") {
  sorted_levels <- eval(substitute(
    df[df$Model == by_model & df$Property == by_prop, ] %>%
      dplyr::group_by(frag) %>%
      dplyr::summarise(res = FUN(Contribution)) %>%
      dplyr::arrange(res),
    list(frag = as.name(frag_col_name))))
  df[, frag_col_name] <- factor(df[, frag_col_name],
                                levels = rspci:::get_col(sorted_levels, 1))
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
  s$ptext <- sapply(rspci:::replace_na(s$ptext), function(i) paste0(rep("*", i), collapse = ""))
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
#' @param x.labels.angle rorate axis x text labels. It works only if \code{flip} was set to TRUE.
#' @return ggplot object which can be further modified before printing
#' @export
#' @examples
#' file_name <- system.file("extdata", "free-wilson_frag_contributions.txt", package = "rspci")
#' df <- load_data(file_name)
#' df <- reorder_data(df, "consensus", "overall")
#' df <- add_signif(df)
#' plot_contrib(df)
#' plot_contrib(df, plot_type = "barplot", flip = FALSE, show_sign_text = "ptext")
plot_contrib <- function(df, frag_name_col = "full_name", contrib_col = "Contribution",
                         plot_type = "boxplot", FUN = median, show_sign_text = NULL,
                         show_sep_lines = TRUE, flip = TRUE, x_labels_angle = 60) {

  # load and detach ggplot2 package if it was not loaded by user
  loaded <- "package:ggplot2" %in% search()
  if (!loaded) {
    if(!require(ggplot2)) {
      return()
    }
  }

  tryCatch({

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
#
#         if (flip) {
#           g <- g + stat_summary(fun.y = FUN, geom = "text", position = position_dodge(width = 0.9), angle = 90,
#                                 aes_string(group = "Property", label = show_sign_text))
#         } else {
#           g <- g + stat_summary(fun.y = FUN, geom = "text", position = position_dodge(width = 0.9),
#                                 aes_string(group = "Property", label = show_sign_text))
#         }
#       } else {
#         if (plot_type == "boxplot") {
#           if ((!flip & length(df$Property) > 1) | (flip & length(df$Property) == 1)) {
#             angle <- 90
#           } else {
#             angle <- 0
#           }
#           g <- g + geom_text(aes_string(group = "Property", label = show_sign_text, y = max(df[, contrib_col])),
#                              position = position_dodge(width = 0.9), angle = angle)

#           if (flip) {
#             g <- g + geom_text(aes_string(group = "Property", label = show_sign_text, y = max(df[, contrib_col])),
#                                position = position_dodge(width = 0.9), angle = 90)
#           } else {
#             g <- g + geom_text(aes_string(group = "Property", label = show_sign_text, y = max(df[, contrib_col])),
#                                position = position_dodge(width = 0.9))
#           }
#         }
#       }
#     }

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

  }, finally = {if (!loaded) detach(package:ggplot2)})

}
