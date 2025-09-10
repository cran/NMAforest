#' Generate Forest Plot with Proportion Contributions for Network Meta-analysis
#' @description
#' This function extends the capabilities of network meta-analysis visualization by
#' generating a forest plot that incorporates direct, indirect, and network meta-analysis (NMA)
#' treatment effects, along with contribution proportions from individual studies and comparison paths.
#' It builds on the statistical foundation provided by the 'netmeta' package and is particularly useful
#' for evaluating the influence of study evidence in complex treatment networks.
#' @details
#' This package includes an implementation of `comparisonStreams()` originally
#' developed by Papakonstantinou et al. (2018) and released with the paper's supplementary material.
#'
#' @param data A data frame in long format with required variables
#'   \code{study} (study identifier) and \code{t} (treatment label), plus outcome
#'   columns: for binary outcomes, \code{r} (events) and \code{N} (sample size);
#'   for continuous outcomes, \code{y} (mean), \code{sd} (standard deviation), and
#'   \code{N} (sample size). Column names are configurable via \code{study},
#'   \code{treat}, \code{event}, \code{N}, \code{mean}, and \code{sd}.
#' @param sm A Character string specifying the summary measure to use (e.g., "OR" for odds ratio, "RR" for risk ratio, "RD" for risk difference, "MD" for mean difference, "SMD" for standardized mean difference, etc.).
#' @param reference Specify the reference treatment for comparisons.
#' @param model Choose "random" or "fixed" effect model.
#' @param comparison A vector of two treatments to compare.
#' @param study Column name identifying the study label.
#' @param treat Column name for treatment assignment.
#' @param event Column name for event counts (for binary outcomes).
#' @param N Column name for total sample size (for binary and continuous outcomes).
#' @param mean Column name for mean values (for continuous outcomes).
#' @param sd Column name for standard deviations (for continuous outcomes).
#' @param study_id Column name used to uniquely identify each study arm (default = "study_id"). If the column is not present in the original dataset, the updated data frame will be returned with this column added.
#' @param study_path Logical. TRUE to include study combination contributions, FALSE to exclude.
#' @return A list containing:
#' \describe{
#'   \item{plot}{A ggplot object of the forest plot.}
#'   \item{output}{A data frame summarizing the results for the specified treatment
#'   comparison, including study-level and overall estimates. Columns contain the
#'   effect sizes with their 95\% confidence intervals, standard errors, study
#'   weights, and contribution proportions, along with labels indicating whether
#'   the row corresponds to a study, direct effect, indirect effect, or the
#'   overall NMA estimate.}
#'   \item{updated_df}{The input data frame after preprocessing, with a unique
#'   study identifier column (from the argument \code{study_id}) added when not
#'   already present. This guarantees consistent referencing of studies in both
#'   the analysis and the output.}
#' }

#' @importFrom netmeta netmeta hatmatrix netsplit
#' @importFrom meta pairwise
#' @importFrom magrittr %>%
#' @importFrom dplyr filter mutate select bind_rows left_join row_number where
#' @importFrom igraph graph_from_data_frame as_ids all_simple_paths set.vertex.attribute graph_from_edgelist set.edge.attribute edge_attr get.shortest.paths delete_edges head_of tail_of E V E<- V<-
#' @importFrom tibble tibble as_tibble
#' @importFrom ggplot2 ggplot aes annotate geom_point guides geom_errorbarh scale_size_continuous geom_segment geom_text geom_vline scale_color_manual scale_y_continuous expand_limits scale_x_continuous labs annotation_custom theme_minimal theme element_text element_line element_blank
#' @importFrom scales rescale
#' @importFrom grid textGrob gpar
#' @importFrom utils head
#' @importFrom stats setNames
#' @importFrom rlist list.append

#' @export
#' @examples
#' # Load example data from the package
#' data(example_data)

#' # Generate the forest plot and return effect estimates and proportion contributions

#' NMAforest(
#'   data = example_data,
#'   sm = "OR",
#'   reference = "x",
#'   model = "random",
#'   comparison = c("x", "y"),
#'   study = "study",
#'   treat = "t",
#'   event = "r",
#'   N = "n",
#'   study_id = "id",
#'   study_path = TRUE
#' )
#' @references
#' Balduzzi, S., Rücker, G., Nikolakopoulou, A., Papakonstantinou, T., Salanti, G., Efthimiou, O., & Schwarzer, G. (2023).
#' netmeta: An R Package for Network Meta-Analysis Using Frequentist Methods. *Journal of Statistical Software*, 106(2), 1–40.
#' https://doi.org/10.18637/jss.v106.i02
#'
#' Csardi, G., & Nepusz, T. (2006). The igraph software package for complex network research. *InterJournal*, Complex Systems, 1695.
#' https://igraph.org
#'
#' Csárdi, G., Nepusz, T., Traag, V., Horvát, Sz., Zanini, F., Noom, D., & Müller, K. (2025).
#' *igraph: Network Analysis and Visualization in R*. https://doi.org/10.5281/zenodo.7682609
#'
#' Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis*. Springer-Verlag New York.
#' https://ggplot2.tidyverse.org
#'
#' Papakonstantinou, T., Nikolakopoulou, A., Rücker, G., Chaimani, A., Schwarzer, G., Egger, M., & Salanti, G. (2018).
#' Estimating the contribution of studies in network meta-analysis: paths, flows and streams. *F1000Research*, 7:610.
#' https://doi.org/10.12688/f1000research.14527.3
#'
#' flow_contribution GitHub repository: https://github.com/esm-ispm-unibe-ch/flow_contribution


# ------------------------------------------------------------------------------
# This file includes code adapted from:
#   streamstatistics2.R in the flow_contribution package
#   https://github.com/esm-ispm-unibe-ch/flow_contribution
#   License: GPL-2
# ------------------------------------------------------------------------------


NMAforest <- function(data,
                      sm,
                      reference,
                      model = "random",
                      comparison,
                      study = "study",
                      treat = "t",
                      event = "r",
                      N = "n",
                      mean = "y",
                      sd = "sd",
                      study_id = "study_id",
                      study_path = TRUE)
{
  
  df <- data
  
  # Generate id column if missing
  
  id_variants <- tolower(names(df))
  has_id <- any(id_variants %in% c("id", "study_id", "studyid", "row_id", "row id", "rowid"))
  
  if (!has_id){
    study_map <- data.frame(
      .study = unique(df[[study]]),
      .study_id = seq_along(unique(df[[study]]))
    )
    names(study_map) <- c(study, study_id)
    df <- left_join(df, study_map, by = study)
  }
  
  # Generate pairwise format
  if (sm %in% c("OR", "RR", "RD", "ASD")) {
    data_long <- eval(substitute(
      pairwise(
        treat = TREAT,
        event = EVENT,
        n = N_,
        studlab = STUDY,
        data = df,
        sm = sm
      ),
      list(
        TREAT = as.name(treat),
        EVENT = as.name(event),
        N_ = as.name(N),
        STUDY = as.name(study)
      )
    ))
  } else if (sm %in% c("MD", "SMD", "ROM")) {
    data_long <- eval(substitute(
      pairwise(
        treat = TREAT,
        mean = MEAN,
        sd = SD,
        n = N_,
        studlab = STUDY,
        data = df,
        sm = sm
      ),
      list(
        TREAT = as.name(treat),
        MEAN = as.name(mean),
        SD = as.name(sd),
        N_ = as.name(N),
        STUDY = as.name(study)
      )
    ))
  }
  
  
  # netmeta and netsplit
  netmeta_result <- netmeta(
    TE = data_long$TE,
    seTE = data_long$seTE,
    treat1 = data_long$treat1,
    treat2 = data_long$treat2,
    studlab = data_long$studlab,
    sm = sm,
    reference.group = reference,
    random = (model == "random"),
    common = (model == "fixed")
  )
  netsplit_result <- netsplit(netmeta_result)
  
  # Build H matrix
  metaNetw <- netmeta_result
  krahn_fn <- tryCatch(
    get("nma.krahn", envir = asNamespace("netmeta"), inherits = FALSE),
    error = function(e) get("nma_krahn", envir = asNamespace("netmeta"), inherits = FALSE)
  )
  
  krahn <- if (model == "fixed") {
    krahn_fn(metaNetw, tau.preset = 0)
  } else {
    krahn_fn(metaNetw, tau.preset = metaNetw$tau)
  }
  
  Hkrahn <- krahn$H
  X.full <- krahn$X.full
  direct <- krahn$direct
  X <- X.full[direct$comparison, , drop = FALSE]
  Vd <- diag(direct$seTE^2, nrow = length(direct$seTE), ncol = length(direct$seTE))
  H <- X.full %*% solve(t(X) %*% solve(Vd) %*% X) %*% t(X) %*% solve(Vd)
  colnames(H) <- rownames(X)
  tau <- netmeta_result$tau
  
  comp1 <- comparison[1]
  comp2 <- comparison[2]
  comp_name1 <- paste(comp2, comp1, sep = ":")
  comp_name2 <- paste(comp1, comp2, sep = ":")
  
  direct_result <- if (model == "random") netsplit_result$direct.random else netsplit_result$direct.fixed
  indirect_result <- if (model == "random") netsplit_result$indirect.random else netsplit_result$indirect.fixed
  TE_matrix <- if (model == "random") netmeta_result$TE.random else netmeta_result$TE.fixed
  lower_matrix <- if (model == "random") netmeta_result$lower.random else netmeta_result$lower.fixed
  upper_matrix <- if (model == "random") netmeta_result$upper.random else netmeta_result$upper.fixed
  
  # Build plot_data
  plot_data <- tibble(Label = character(), EffectSize = numeric(), LowerCI = numeric(), UpperCI = numeric(), Type = character())
  
  d <- direct_result  %>%
    filter((comparison == comp_name1) | (comparison == comp_name2))
  
  plot_data <- bind_rows(plot_data, tibble(
    Label = "Overall Direct Effect",
    EffectSize = ifelse(d$comparison == comp_name1, d$TE, -d$TE),
    LowerCI = ifelse(d$comparison == comp_name1, d$lower, -d$upper),
    UpperCI = ifelse(d$comparison == comp_name1, d$upper, -d$lower),
    Type = "Overall"
  ))
  
  # Study-level Direct Effect
  direct_rows <- data_long %>%
    filter((treat1 == comp1 & treat2 == comp2) | (treat1 == comp2 & treat2 == comp1))
  
  for (i in seq_len(nrow(direct_rows))) {
    row <- direct_rows[i, ]
    te <- if (row$treat1 == comp1 & row$treat2 == comp2) -row$TE else row$TE
    
    if (model == "random") {
      ci_low <- te - 1.96 * sqrt(row$seTE^2 + tau^2)
      ci_up  <- te + 1.96 * sqrt(row$seTE^2 + tau^2)
    } else {
      ci_low <- te - 1.96 * sqrt(row$seTE^2)
      ci_up  <- te + 1.96 * sqrt(row$seTE^2)
    }
    
    plot_data <- bind_rows(plot_data, tibble(
      Label = paste0("Study ", row[[study_id]]),
      EffectSize = te,
      LowerCI = ci_low,
      UpperCI = ci_up,
      Type = "Direct"
    ))
  }
  
  id <- indirect_result %>%
    filter((comparison == comp_name1) | (comparison == comp_name2))
  
  plot_data <- bind_rows(plot_data, tibble(
    Label = "Overall Indirect Effect",
    EffectSize = ifelse(id$comparison == comp_name1, id$TE, -id$TE),
    LowerCI = ifelse(id$comparison == comp_name1, id$lower, -id$upper),
    UpperCI = ifelse(id$comparison == comp_name1, id$upper, -id$lower),
    Type = "Overall"
  ))
  
  # Indirect Paths
  max_combos_per_path <- 30
  edges <- data_long %>% select(treat1, treat2)
  g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  all_paths <- igraph::all_simple_paths(g, from = comp1, to = comp2, cutoff = 4)
  indirect_paths <- all_paths[sapply(all_paths, length) > 2]
  
  get_theta <- function(from, to, effect_table) {
    comp_forward <- paste(from, to, sep = ":")
    comp_reverse <- paste(to, from, sep = ":")
    row <- effect_table %>% filter(comparison %in% c(comp_forward, comp_reverse))
    if (nrow(row) == 0) return(NA)
    
    if (row$comparison == comp_forward) {
      return(-row$TE)
    } else {
      return(row$TE)
    }
  }
  
  netsplit_direct <- netsplit_result$direct.random
  
  indirect_path_info <- list()
  for (p in indirect_paths) {
    path_nodes <- igraph::as_ids(p)
    path_label <- paste("Path", paste(path_nodes, collapse = "\u2192 "))
    
    # Compute indirect effect and record segment signs
    effect <- 0
    segment_list <- list()
    valid_effect <- TRUE
    
    for (i in seq_len(length(path_nodes) - 1)) {
      from <- path_nodes[i]
      to <- path_nodes[i + 1]
      comp_forward <- paste(from, to, sep = ":")
      comp_reverse <- paste(to, from, sep = ":")
      row <- netsplit_direct %>% filter(comparison %in% c(comp_forward, comp_reverse))
      
      if (nrow(row) == 0) {
        effect <- NA
        valid_effect <- FALSE
        break
      }
      
      if (row$comparison == comp_forward) {
        theta <- -row$TE  # from - to
        segment_sign <- +1
      } else {
        theta <- row$TE  # to - from
        segment_sign <- -1
      }
      
      effect <- effect + theta
      segment_list[[i]] <- list(from = from, to = to, sign = segment_sign)
    }
    
    # CI Calculation for this indirect path
    get_variance <- function(from, to) {
      comp1 <- paste(from, to, sep = ":")
      comp2 <- paste(to, from, sep = ":")
      row <- netsplit_result$direct.random %>% filter(comparison %in% c(comp1, comp2))
      if (nrow(row) == 0) return(NA)
      if (model == "random") {
        return((row$seTE)^2 + tau^2)
      } else {
        return((row$seTE)^2)
      }
    }
    
    # Segment and variance info
    segments <- list()
    segment_vars <- c()
    segment_ses <- c()
    
    for (i in seq_len(length(path_nodes) - 1)) {
      from <- path_nodes[i]
      to <- path_nodes[i + 1]
      var <- get_variance(from, to)
      if (is.na(var)) {
        total_var <- NA
        break
      }
      segments[[i]] <- c(from, to)
      segment_vars[i] <- var
      segment_ses[i] <- sqrt(var)
    }
    
    # Total variance: sum of all variances
    total_var <- sum(segment_vars)
    
    # Covariance terms with fixed sign logic
    rho <- 0.5
    for (i in seq_len(length(path_nodes) - 2)) {
      for (j in (i + 1):(length(path_nodes) - 1)) {
        a1 <- path_nodes[i]
        a2 <- path_nodes[i + 1]
        b1 <- path_nodes[j]
        b2 <- path_nodes[j + 1]
        
        ids_a <- data_long[[study_id]][
          (data_long$treat1 == a1 & data_long$treat2 == a2) |
            (data_long$treat1 == a2 & data_long$treat2 == a1)
        ]
        ids_b <- data_long[[study_id]][
          (data_long$treat1 == b1 & data_long$treat2 == b2) |
            (data_long$treat1 == b2 & data_long$treat2 == b1)
        ]
        
        shared_ids <- intersect(ids_a, ids_b)
        if (length(shared_ids) > 0) {
          se_a <- sqrt(get_variance(a1, a2))
          se_b <- sqrt(get_variance(b1, b2))
          
          if (i == 1 || j == 1) {
            sign_cov <- -1  # interaction of positive and negative
          } else {
            sign_cov <- 1  # both negative
          }
          
          cov_term <- sign_cov * rho * se_a * se_b
          total_var <- total_var + 2 * cov_term
        }
      }
    }
    
    # Confidence interval
    if (!is.na(total_var)) {
      se_indirect <- sqrt(total_var)
      ci_low <- effect - 1.96 * se_indirect
      ci_up  <- effect + 1.96 * se_indirect
    } else {
      ci_low <- NA
      ci_up <- NA
    }
    
    plot_data <- bind_rows(plot_data, tibble(
      Label = path_label,
      EffectSize = effect,
      LowerCI = ci_low,
      UpperCI = ci_up,
      Type = "Indirect"
    ))
    
    indirect_path_info[[path_label]] <- path_nodes
    
    
    # Add Study ID combinations for path
    study_sets <- list()
    valid <- TRUE
    for (i in seq_along(path_nodes)[-length(path_nodes)]) {
      t1 <- path_nodes[i]
      t2 <- path_nodes[i + 1]
      rows <- data.frame(
        id = data_long[[study_id]][
          (data_long$treat1 == t1 & data_long$treat2 == t2) |
            (data_long$treat1 == t2 & data_long$treat2 == t1)
        ]
      )
      if (nrow(rows) == 0) {
        valid <- FALSE
        break
      }
      study_sets[[i]] <- rows
    }
    if (!valid) next
    combos <- study_sets[[1]]
    for (i in 2:length(study_sets)) {
      combos <- merge(combos, study_sets[[i]], by = NULL)
    }
    if (nrow(combos) > max_combos_per_path) {
      combos <- combos[1:max_combos_per_path, ]
    }
    
    if (study_path) {
      for (j in seq_len(nrow(combos))) {
        ids <- unlist(combos[j, grep("^id", names(combos))])
        label <- paste0("Study ", paste(ids, collapse = "-"))
        plot_data <- bind_rows(plot_data, tibble(
          Label = label,
          EffectSize = NA, LowerCI = NA, UpperCI = NA, Type = "Blank"
        ))
      }
    }
  }
  
  # NMA Effect
  if (comp1 %in% rownames(TE_matrix) && comp2 %in% colnames(TE_matrix)) {
    plot_data <- bind_rows(plot_data, tibble(
      Label = "Overall NMA Effect",
      EffectSize = -TE_matrix[comp1, comp2],
      LowerCI = -lower_matrix[comp1, comp2],
      UpperCI = -upper_matrix[comp1, comp2],
      Type = "Overall"
    ))
  }
  
  plot_data <- plot_data %>%
    mutate(OrderIndex = nrow(.) - row_number() + 1)
  
  tau <- netmeta_result$tau
  
  plot_data <- plot_data %>%
    mutate(
      raw_se = (UpperCI - LowerCI) / (2 * 1.96),
      weight = ifelse(!is.na(raw_se), 1 / (raw_se^2), NA)
    )
  
  # Get flows
  comparisonStreams = function(hatmatrix, comparison){
    
    directs <- colnames(hatmatrix)
    
    hatMatrix <- hatmatrix
    
    rownames(hatMatrix) <- rownames(hatmatrix)
    
    split <- function (dir) {strsplit(dir,":")}
    
    dims <- dim(hatMatrix)
    
    #rows of comparison matrix
    comparisons <- unlist(lapply(rownames(hatMatrix),unlist))
    
    comparisonToEdge <- function (comp) unlist (split(comp))
    
    directlist <- unlist(lapply(lapply(directs,split),unlist))
    # print(c("dir",directs))
    
    edgeList <- matrix( directlist, ncol = 2, byrow = TRUE)
    # print(c("Edgelist"))
    # print(edgeList)
    
    g <- graph_from_edgelist(edgeList , directed=FALSE)
    g <- set.vertex.attribute(g,'label',value = V(g))
    g <- set.edge.attribute(g,'label',value = E(g))
    #print(V(g)$label)
    #print(V(g)$name)
    #print(E(g))
    
    setWeights <- function (g,comparison,conMat) {
      set.edge.attribute(g,"weight",value=rep(0,dims[2]))
    }
    
    
    getFlow <- function(g,edge) {return(E(g)[edge]$flow)}
    
    sv <- function (comparison) {split(comparison)[[1]][1][1]}
    
    tv <- function (comparison) {split(comparison)[[1]][2][1]}
    
    initRowGraph <- function(comparison) {
      dedgeList <- lapply(1:length(directs),function(comp) {
        if(hatMatrix[comparison,comp]>0){
          # print(c("not switched",directs[comp],hatMatrix[comparison,comp]))
          return (c(sv(directs[comp]),tv(directs[comp])))
        }else{
          # print(c("switched",directs[comp],hatMatrix[comparison,comp]))
          return (c(tv(directs[comp]),sv(directs[comp])))
        }
      })
      dedgeList <- matrix( unlist(dedgeList), ncol = 2, byrow = TRUE)
      # gg <- setFlow(g,comparison)
      # E(gg)$weight <- rep(0,dims[2])
      # return(gg)
      flows<-abs(hatMatrix[comparison,])
      dg <- graph_from_edgelist(dedgeList , directed = TRUE)
      E(dg)[]$weight <- rep(0,dims[2])
      E(dg)[]$flow <- abs(hatMatrix[comparison,])
      V(dg)[]$label <- V(dg)[]$name
      # E(dg)[]$label <- E(dg)[]$flow
      dg <- set.edge.attribute(dg,'label',value = E(dg))
      # print(c("isdirected",is.directed(dg)))
      return(dg)
    }
    
    contribution = rep(0,dims[2])
    streams = list()
    names(contribution) <- c(1:dims[2])
    
    reducePath <- function (g,comparison,spl) {
      pl <- length(spl[[1]])
      splE <- lapply(spl[[1]], function(e){
        return (E(g)[e[]])
      })
      flow <- min(unlist(lapply(splE, function(e){
        return(e$flow[])
      })))
      path = toString(names(unlist(lapply(spl,function(e){c(head_of(g,e),tail_of(g,e))}))))
      streams <<- list.append(streams,data.frame(comp=comparison,length=floor(length(splE)),stream=path,flow=flow))
      # print(c("to shortest path einai :",spl))
      gg <- Reduce(function(g, e){
        elabel <- e$label
        # print(c("pame plevra:",e,"dld",e$label))
        pfl <- e$flow[]
        g <- set.edge.attribute(g,"flow",e, pfl-flow)
        # print(c("h e",e,"einai pragmatika h ",elabel))
        cw <-  e$weight[] + (flow[1]/pl)
        # print(c("flow",flow,"eweight",e$weight[]))
        contribution[elabel] <<- cw
        return(set.edge.attribute(g,"weight",e, cw))
      },splE, g)
      # print(c("graph before deleting edges", E(gg)$label))
      emptyEdges <- Reduce(function(removedEdges, e){
        e <- E(gg)[e[]]
        if(e$flow[[1]][[1]]==0){
          removedEdges <- c(removedEdges, e)
        }
        return(removedEdges)
      },splE, c())
      # print(c("edges to be removed",emptyEdges))
      return(delete_edges(gg, emptyEdges))
      # print(c("graph after deleting edges", E(gg)$label))
    }
    reduceGraph <- function (g,comparison) {
      getshortest <- function (g,compariston) {
        floweights = lapply(edge_attr(g,"flow",E(g)), function(f){return(abs(2-f))})
        spths = suppressWarnings(
          get.shortest.paths(g,sv(comparison),tv(comparison),mode="out",output="epath",weights=floweights)
        )
        return(spths$epath)
      }
      # while(edge_connectivity(g,sv(comparison),tv(comparison))>0){
      spath <- getshortest(g,comparison)
      while(length(unlist(spath))>0){
        g <- reducePath(g,comparison,spath)
        spath <- getshortest(g,comparison)
      }
      # print("teleiwse")
      return(g)
    }
    
    # ptm <- proc.time()
    # gg <- reduceGraph (initRowGraph(comparison), comparison)
    reduceGraph (initRowGraph(comparison), comparison)
    # executionTime <- proc.time() - ptm
    # print(c("execution time",executionTime))
    
    names(contribution) <- directs
    contribution <- 100 * contribution
    
    return(list( streams=streams
                 ,contribution=contribution
                 
    ))
  }
  
  # Build proportion bars
  sorted_comp1 <- paste(comp1, comp2, sep = ":")
  sorted_comp2 <- paste(comp2, comp1, sep = ":")
  rowname_match <- rownames(H)[rownames(H) %in% c(sorted_comp1, sorted_comp2)]
  colname_match <- colnames(H)[colnames(H) %in% c(sorted_comp1, sorted_comp2)]
  overall_direct_prop <- round(abs(H[rowname_match, colname_match]), 3)
  w <- 1 / (netmeta_result$seTE.adj.common^2 + tau^2)
  study_rows <- which((data_long$treat1 == comp1 & data_long$treat2 == comp2) |
                        (data_long$treat1 == comp2 & data_long$treat2 == comp1))
  # Calculate direct bar proportions safely
  if (length(study_rows) > 0) {
    study_weights <- w[study_rows]
    study_props <- study_weights / sum(study_weights) * overall_direct_prop
    study_ids <- data_long[[study_id]][study_rows]
    cumulative_start <- cumsum(c(0, head(study_props, -1))) + (max(plot_data$UpperCI, na.rm = TRUE) + 0.5)
    cumulative_end <- cumulative_start + study_props
    
    direct_bar_data <- tibble(
      xstart = cumulative_start,
      xend = cumulative_end,
      y = plot_data$OrderIndex[match(paste0("Study ", study_ids), plot_data$Label)],
      Type = "Direct",
      Proportion = round(study_props, 3)
    )
    
    overall_bar <- tibble(
      xstart = min(cumulative_start),
      xend = max(cumulative_end),
      y = plot_data$OrderIndex[plot_data$Label == "Overall Direct Effect"],
      Type = "Overall",
      Proportion = overall_direct_prop
    )
  } else {
    # No direct comparisons
    cumulative_start <- numeric(0)
    cumulative_end <- numeric(0)
    direct_bar_data <- tibble()
    overall_bar <- tibble()
  }
  
  
  # Safely define overall_direct_prop, fallback to NA if not present
  overall_direct_prop <- tryCatch({
    if (exists("overall_direct_prop") && length(overall_direct_prop) > 0) {
      overall_direct_prop
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)
  
  # Determine xstart for indirect bar
  xstart_indirect <- if (length(cumulative_end) > 0 && all(is.finite(cumulative_end))) {
    max(cumulative_end, na.rm = TRUE)
  } else {
    max(plot_data$UpperCI, na.rm = TRUE) + 0.2
  }
  
  # Compute indirect proportion and create bar
  overall_indirect_prop <- if (is.na(overall_direct_prop)) {
    1
  } else {
    round(1 - overall_direct_prop, 3)
  }
  
  overall_indirect_bar <- tibble(
    xstart = xstart_indirect,
    xend = xstart_indirect + overall_indirect_prop,
    y = plot_data$OrderIndex[plot_data$Label == "Overall Indirect Effect"],
    Type = "Overall",
    Proportion = overall_indirect_prop
  )
  
  streams_out <- comparisonStreams(H, ifelse(any(comp_name1 %in% rownames(H)), comp_name1, comp_name2))$streams
  
  cumulative_start_indirect <- if (exists("cumulative_end") && length(cumulative_end) > 0 && all(is.finite(cumulative_end))) {
    max(cumulative_end)
  } else {
    max(plot_data$UpperCI, na.rm = TRUE) + 0.2
  }
  
  #proportion bar
  path_node_map <- lapply(indirect_path_info, function(nodes) sort(unique(nodes)))
  grouped_paths <- split(names(path_node_map), sapply(path_node_map, paste, collapse = "-"))
  path_divisor <- sapply(grouped_paths, length)
  path_group_lookup <- unlist(lapply(names(grouped_paths), function(k) {
    setNames(rep(length(grouped_paths[[k]]), length(grouped_paths[[k]])), grouped_paths[[k]])
  }))
  
  indirect_bar_data <- tibble()
  
  for (path_label in names(indirect_path_info)) {
    path_nodes <- indirect_path_info[[path_label]]
    path_length <- length(path_nodes) - 1
    path_label_clean <- gsub("\\s*\u2192\\s*", "\u2192", path_label)
    plot_labels_clean <- gsub("\\s*\u2192\\s*", "\u2192", plot_data$Label)
    y_value <- plot_data$OrderIndex[match(path_label_clean, plot_labels_clean)]
    
    matched <- FALSE
    
    for (stream in streams_out) {
      stream_nodes <- unlist(strsplit(as.character(stream$stream), ","))
      stream_nodes <- trimws(stream_nodes)  # remove spaces
      stream_nodes <- stream_nodes[stream_nodes != ""]  # remove empty
      
      stream_nodes_clean <- unique(stream_nodes)
      path_nodes_clean <- unique(path_nodes)
      
      if (length(path_nodes_clean) == length(stream_nodes_clean)) {
        if (all(path_nodes_clean %in% stream_nodes_clean) && all(stream_nodes_clean %in% path_nodes_clean)) {
          
          flow_val <- stream$flow/ path_group_lookup[[path_label]]
          indirect_bar_data <- bind_rows(indirect_bar_data, tibble(
            xstart = cumulative_start_indirect,
            xend = cumulative_start_indirect + flow_val,
            y = y_value,
            Type = "Indirect",
            Proportion = round(flow_val, 3)
          ))
          cumulative_start_indirect <- cumulative_start_indirect + flow_val
          matched <- TRUE
          break  # matched, go to next path
        }
      }
    }
    
    if (!matched) {
      warning(paste0("No stream matched for path: ", path_label))
    }
    
    # Build sub-bars for study combinations within each indirect path
    sub_bar_data <- tibble()
    
    for (path_label in names(indirect_path_info)) {
      path_nodes <- indirect_path_info[[path_label]]
      path_label_clean <- gsub("\\s*\u2192\\s*", "\u2192", path_label)
      plot_labels_clean <- gsub("\\s*\u2192\\s*", "\u2192", plot_data$Label)
      y_path_value <- plot_data$OrderIndex[match(path_label_clean, plot_labels_clean)]
      
      indirect_bar_row <- indirect_bar_data %>% filter(y == y_path_value)
      if (nrow(indirect_bar_row) == 0) next
      
      indirect_start <- indirect_bar_row$xstart
      indirect_flow <- indirect_bar_row$xend - indirect_bar_row$xstart
      
      # Rebuild study_sets
      study_sets <- list()
      valid <- TRUE
      for (i in seq_along(path_nodes)[-length(path_nodes)]) {
        t1 <- path_nodes[i]
        t2 <- path_nodes[i + 1]
        comp_label <- paste(pmin(t1, t2), pmax(t1, t2), sep = ":")
        
        rows <- data.frame(
          id = data_long[[study_id]][
            (data_long$treat1 == t1 & data_long$treat2 == t2) |
              (data_long$treat1 == t2 & data_long$treat2 == t1)
          ],
          comp = comp_label
        )
        
        if (nrow(rows) == 0) {
          valid <- FALSE
          break
        }
        study_sets[[i]] <- rows
      }
      if (!valid) next
      
      combos <- study_sets[[1]]
      for (i in 2:length(study_sets)) {
        combos <- merge(combos, study_sets[[i]], by = NULL)
      }
      if (nrow(combos) > max_combos_per_path) {
        combos <- combos[1:max_combos_per_path, ]
      }
      
      cumulative_substart <- indirect_start
      
      for (j in seq_len(nrow(combos))) {
        ids <- unlist(combos[j, grep("^id", names(combos))])
        comps <- unlist(combos[j, grep("^comp", names(combos))])
        if (length(ids) != length(comps)) next
        
        partials <- numeric(length(ids))
        for (k in seq_along(ids)) {
          id_k <- ids[k]
          comp_k <- comps[k]
          t1 <- strsplit(comp_k, ":")[[1]][1]
          t2 <- strsplit(comp_k, ":")[[1]][2]
          
          rows_k <- which((data_long$treat1 == t1 & data_long$treat2 == t2) |
                            (data_long$treat1 == t2 & data_long$treat2 == t1))
          if (length(rows_k) == 0) next
          
          weights_k <- w[rows_k]
          this_weight <- w[data_long[[study_id]] == id_k &
                             ((data_long$treat1 == t1 & data_long$treat2 == t2) |
                                (data_long$treat1 == t2 & data_long$treat2 == t1))]
          
          if (length(this_weight) != 1 || sum(weights_k) == 0) {
            partials[k] <- NA
          } else {
            partials[k] <- this_weight / sum(weights_k)
          }
        }
        
        final_prop <- prod(partials, na.rm = TRUE) * indirect_flow
        if (!is.na(final_prop) && final_prop > 0) {
          label_combination <- paste0("Study ", paste(ids, collapse = "-"))
          
          label_matches <- which(plot_data$Label == label_combination)
          y_options <- plot_data$OrderIndex[label_matches]
          y_match <- y_options[y_options < y_path_value]
          
          y <- if (length(y_match) > 0) max(y_match) else NA
          if (!is.na(y)) {
            sub_bar_data <- bind_rows(sub_bar_data, tibble(
              xstart = cumulative_substart,
              xend = cumulative_substart + final_prop,
              y = y,
              Type = "Blank",
              Proportion = round(final_prop, 3)
            ))
            cumulative_substart <- cumulative_substart + final_prop
          }
        }
      }
    }
  }
  
  xstart_total <- if (length(cumulative_start) > 0) {
    min(cumulative_start)
  } else if (length(cumulative_end) > 0) {
    max(cumulative_end)
  } else {
    max(plot_data$UpperCI, na.rm = TRUE) + 0.2
  }
  
  xend_total <- if (nrow(indirect_bar_data) > 0 && all(is.finite(indirect_bar_data$xend))) {
    max(indirect_bar_data$xend)
  } else {
    xstart_total + 1
  }
  
  
  total_bar <- tibble(
    xstart = xstart_total,
    xend = xend_total,
    y = plot_data$OrderIndex[plot_data$Label == "Overall NMA Effect"],
    Type = "Overall",
    Proportion = 1
  )
  
  bar_data <- bind_rows(overall_bar, direct_bar_data, overall_indirect_bar, indirect_bar_data, sub_bar_data, total_bar)
  
  bar_data <- bind_rows(
    if (nrow(overall_bar) > 0) overall_bar else NULL,
    if (nrow(direct_bar_data) > 0) direct_bar_data else NULL,
    overall_indirect_bar,
    indirect_bar_data,
    sub_bar_data,
    total_bar
  )
  
  right_edge <- max(bar_data$xend, na.rm = TRUE)
  bar_data <- bar_data %>%
    mutate(x_label = right_edge)
  
  # First compute the max x position with grid (effect size only)
  x_grid_max <- max(plot_data$UpperCI, na.rm = TRUE) + 0.2
  x_text_pos <- max(bar_data$xend, na.rm = TRUE) + 0.05
  vline_start <- min(bar_data$xstart, na.rm = TRUE)
  vline_end_all <- max(bar_data$xend, na.rm = TRUE)
  vline_end_direct <- bar_data$xend[bar_data$y == plot_data$OrderIndex[plot_data$Label == "Overall Direct Effect"]]
  
  if (all(is.na(plot_data$LowerCI)) || all(is.na(plot_data$UpperCI))) {
    effect_center <- 0
  } else {
    effect_center <- mean(c(
      min(plot_data$LowerCI, na.rm = TRUE),
      max(plot_data$UpperCI, na.rm = TRUE)
    ))
  }
  
  if (!"xstart" %in% names(bar_data) || !"xend" %in% names(bar_data)) {
    proportion_center <- xstart_total + 0.5
  } else {
    valid_xstart <- bar_data$xstart[is.finite(bar_data$xstart)]
    valid_xend <- bar_data$xend[is.finite(bar_data$xend)]
    
    if (length(valid_xstart) == 0 || length(valid_xend) == 0) {
      proportion_center <- xstart_total + 0.5
    } else {
      proportion_center <- mean(c(min(valid_xstart), max(valid_xend)))
    }
  }
  
  
  # Estimate the top y-value
  top_y <- max(plot_data$OrderIndex, na.rm = TRUE) + 1.5
  
  bar_one <- bar_data %>%
    dplyr::arrange(y, dplyr::desc(Type != "Blank")) %>%  
    dplyr::distinct(y, .keep_all = TRUE) %>%
    dplyr::select(
      y,
      bar_xstart = xstart,
      bar_xend   = xend,
      bar_prop   = Proportion,
      bar_type   = Type
    )
  
  combined <- plot_data %>%
    dplyr::left_join(bar_one, by = c("OrderIndex" = "y")) %>%
    dplyr::mutate(
      size_scaled = dplyr::if_else(
        !is.na(EffectSize) & !is.na(LowerCI) & !is.na(UpperCI) & !is.na(bar_prop),
        scales::rescale(bar_prop, to = c(4, 9)),
        3
      )
    )
  
  # Construct ggplot
  plot_object <- ggplot(plot_data, aes(y = OrderIndex)) +
    annotate("rect",
             xmin = x_grid_max, xmax = Inf,
             ymin = -Inf, ymax = Inf,
             fill = "white", alpha = 1) +
    
    geom_point(data = subset(combined, !is.na(EffectSize) & Type == "Overall"),
               aes(x = EffectSize, color = Type, size = bar_prop), shape = 19) +
    geom_point(data = subset(combined, !is.na(EffectSize) & Type == "Direct"),
               aes(x = EffectSize, color = Type, size = bar_prop), shape = 18) +
    geom_point(data = subset(combined, !is.na(EffectSize) & Type == "Indirect"),
               aes(x = EffectSize, color = Type, size = bar_prop), shape = 15) +
    scale_size_continuous(range = c(2, 8.5)) +
    
    guides(size = "none") +
    
    geom_errorbarh(data = subset(plot_data, !is.na(LowerCI)),
                   aes(xmin = LowerCI, xmax = UpperCI, color = Type),
                   height = 0.3, size = 0.8) +
    
    geom_segment(data = bar_data,
                 aes(x = xstart, xend = xend, y = y, yend = y),
                 inherit.aes = FALSE, size = 1.8, color = "black") +
    
    geom_text(data = bar_data,
              aes(x = x_text_pos, y = y, label = Proportion),
              inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = 3.5, color = "black") +
    
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_vline(xintercept = vline_start, linetype = "dashed", color = "black") +
    geom_vline(xintercept = vline_end_direct, linetype = "dashed", color = "black") +
    geom_vline(xintercept = vline_end_all, linetype = "dashed", color = "black") +
    
    scale_color_manual(values = c("Direct" = "blue", "Indirect" = "darkorange", "Overall" = "#009E73")) +
    scale_y_continuous(breaks = plot_data$OrderIndex, labels = plot_data$Label) +
    expand_limits(y = top_y) +
    scale_x_continuous(
      breaks = seq(floor(min(plot_data$LowerCI, na.rm = TRUE)),
                   floor(min(bar_data$xstart, na.rm = TRUE)),
                   by = 1)
    ) +
    
    labs(
      title = paste(
        "Forest Plot for Comparison", comp1, "&", comp2,
        "-", ifelse(model == "random", "Random Effects Model", "Fixed Effects Model")
      ),
      x = sm,
      y = NULL,
      color = "Type"
    ) +
    
    annotation_custom(
      grob = grid::textGrob("Effect Size", gp = grid::gpar(fontface = "bold", fontsize = 12)),
      xmin = effect_center,
      xmax = effect_center,
      ymin = top_y,
      ymax = top_y
    ) +
    
    annotation_custom(
      grob = grid::textGrob("Proportion", gp = grid::gpar(fontface = "bold", fontsize = 12)),
      xmin = proportion_center,
      xmax = proportion_center,
      ymin = top_y,
      ymax = top_y
    ) +
    
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "top",
      panel.grid.major.x = element_line(color = "gray80"),
      panel.grid.minor = element_blank()
    )
  
  result <- combined %>%
    dplyr::select(
      Label, EffectSize, LowerCI, UpperCI, SE = raw_se, weight, 
      proportion = bar_prop
    )
  
  result <- result %>%
    dplyr::mutate(
      dplyr::across(
        where(is.numeric),        
        ~ round(., 3)             
      )
    )
  result <- as.data.frame(result)
  result[is.na(result)] <- ""
  
  return(list(
    plot = plot_object,
    output = result,
    updated_df = if (!has_id) df else NULL
  ))
  
}

utils::globalVariables(c(
  "treat1", "treat2", ".", "UpperCI", "LowerCI", "raw_se", "weight", "Label",
  "OrderIndex", "EffectSize", "Type", "size_scaled", "xstart", "xend", "Proportion", "bar_prop"
))

