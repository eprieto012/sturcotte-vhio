# --- Configure Python environment and verify MAGIC availability ---
library(reticulate)
reticulate::use_python(
  "/mnt/bioinfnas/immuno/eprieto/Rmagic/miniforge3/envs/magic_env/bin/python",
  required = TRUE
)
py_config()
py_module_available("magic")
# --- Load R libraries ---
library(data.table)
library(djvdj)
library(dplyr)
library(ggalluvial)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggrepel)
library(ggridges)
library(ggVennDiagram)
library(keras3)
library(pals)
library(patchwork)
library(pdfCluster)
library(ProjecTILs)
library(purrr)
library(RColorBrewer)
library(scGate)
reticulate::import("magic")
library(Rmagic)
library(readr)
library(scDblFinder)
library(Seurat)
library(SoupX)
library(stringr)
library(tensorflow)
library(tidyr)
library(viridis)
# --- Phenotype colours ---
phenotype_colors <- c(
  "CD8.NaiveLike" = "gray80",
  "CD4.NaiveLike" = "gray80",
  "CD8.MAIT" = "#ADFF2F",
  "CD8.CM" = "#f4c40f",
  "CD8.EM" = "#FF69B4",
  "CD8.TEMRA" = "#8A2BE2",
  "CD8.TEX" = "black",
  "CD8.TPEX" = "#fe9b60",
  "CD4.CTL_EOMES" = "#32CD32",
  "CD4.CTL_Exh" = "#4A4999",
  "CD4.CTL_GNLY" = "#7FFFD4",
  "CD4.Tfh" = "#1E90FF",
  "CD4.Th17" = "#556B2F",
  "CD4.Treg" = "#d8443c"
)
# --- Set random seeds in R and Python for reproducibility ---
set.seed(123)
py_run_string("
import numpy as np
import random
np.random.seed(123)
random.seed(123)
")
# --- Configure plotting style and visualization settings ---
theme_ic <- theme_bw(base_size = 11) +
  theme(
    aspect.ratio = 0.5,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(color = "black", face = "plain"),
    plot.title = element_text(hjust = 0),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black")
  )
# --- Define output directories and file paths ---
output_path <- file.path(datapath, sample, folder_data)
save_paths <- list(
  plots = file.path(output_path, "plots"),
  tcr = file.path(output_path, "tcr"),
  markers = file.path(output_path, "markers"),
  seurat = file.path(output_path, "seurat")
)
invisible(
  lapply(
    save_paths,
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  )
)
# --- Configure plot export format and resolution ---
plot_format <- "png"
plot_dpi <- 600
plots_list <- list()
pdfcluster_plots <- list()
tables_list <- list(
  qc = list(),
  frequencies = list(),
  tcr = list(),
  markers = list()
)
# --- Output helpers ---
count_plot_panels <- function(plot) {
  if (inherits(plot, "patchwork")) {
    return(
      tryCatch(
        max(1L, length(plot$patches$plots) + 1L),
        error = function(e) 1L
      )
    )
  }
  if (inherits(plot, "ggplot")) {
    return(
      tryCatch(
        max(
          1L,
          length(
            unique(
              ggplot2::ggplot_build(plot)$layout$layout$PANEL
            )
          )
        ),
        error = function(e) 1L
      )
    )
  }
  1L
}
get_plot_dimensions <- function(
    plot_name,
    plot
) {
  n_panels <- count_plot_panels(plot)
  # --- Larger dimensions for LS11 / IC11 prioritized TCR plots ---
  if (grepl("^(LS11|IC11)_", plot_name)) {
    ncol <- 3L
    nrow <- ceiling(n_panels / ncol)
    
    return(
      c(
        width = 30,
        height = max(10, 6 * nrow)
      )
    )
  }
  if (
    grepl(
      "_All_conditions_genes_magic.*_umap$",
      plot_name
    )
  ) {
    ncol <- max(
      1L,
      ceiling(n_panels / 2)
    )
    return(
      c(
        width = max(22, 5 * ncol),
        height = 9
      )
    )
  }
  if (grepl("_rank$", plot_name)) {
    ncol <- 2L
    nrow <- ceiling(n_panels / ncol)
    return(
      c(
        width = 16,
        height = max(10, 2.5 * nrow)
      )
    )
  }
  if (grepl("_heatmap$", plot_name)) {
    return(
      c(
        width = 16,
        height = 14
      )
    )
  }
  if (grepl("_umap$", plot_name)) {
    return(
      c(
        width = 8,
        height = 8
      )
    )
  }
  if (grepl("_Volcano$", plot_name)) {
    return(
      c(
        width = 8,
        height = 8
      )
    )
  }
  if (grepl("global_TCR_Venn", plot_name)) {
    return(
      c(
        width = 10,
        height = 8
      )
    )
  }
  if (grepl("_top(70|100)$", plot_name)) {
    return(
      c(
        width = 12,
        height = 8
      )
    )
  }
  if (identical(plot_name, "IC07_celltype")) {
    return(
      c(
        width = 14,
        height = 10
      )
    )
  }
  if (grepl("only_vs_non_only", plot_name)) {
    ncol <- min(10L, n_panels)
    nrow <- ceiling(n_panels / ncol)
    return(
      c(
        width = min(24, max(12, 2.4 * ncol)),
        height = min(24, max(7, 3 * nrow + 2))
      )
    )
  }
  if (
    grepl(
      "frequency_FC_prioritization_scatter",
      plot_name
    )
  ) {
    ncol <- min(7L, n_panels)
    nrow <- ceiling(n_panels / ncol)
    return(
      c(
        width = min(24, max(10, 3 * ncol)),
        height = min(24, max(7, 3 * nrow + 2))
      )
    )
  }
  if (grepl("frequency_FC", plot_name)) {
    ncol <- min(10L, n_panels)
    nrow <- ceiling(n_panels / ncol)
    return(
      c(
        width = min(24, max(10, 2.4 * ncol)),
        height = min(24, max(7, 3 * nrow + 2))
      )
    )
  }
  if (startsWith(plot_name, "IC")) {
    ncol <- min(
      4L,
      max(
        1L,
        ceiling(sqrt(n_panels))
      )
    )
    nrow <- ceiling(n_panels / ncol)
    return(
      c(
        width = min(20, max(8, 3.5 * ncol)),
        height = min(20, max(6, 3.2 * nrow + 1))
      )
    )
  }
  ncol <- min(
    4L,
    max(
      1L,
      ceiling(sqrt(n_panels))
    )
  )
  nrow <- ceiling(n_panels / ncol)
  c(
    width = min(20, max(8, 3.5 * ncol)),
    height = min(20, max(6, 3.2 * nrow + 1))
  )
}
# --- Custom R functions ---
save_plot_output <- function(plot_name, plot) {
  dir.create(save_paths$plots, recursive = TRUE, showWarnings = FALSE)
  dimensions <- get_plot_dimensions(plot_name = plot_name, plot = plot)
  ggsave_args <- list(
    filename = file.path(save_paths$plots, paste0(plot_name, ".", plot_format)),
    plot = plot,
    device = plot_format,
    width = dimensions[["width"]],
    height = dimensions[["height"]],
    units = "in",
    dpi = plot_dpi,
    limitsize = FALSE,
    bg = "white"
  )
  if (identical(plot_format, "tiff")) ggsave_args$compression <- "lzw"
  do.call(ggplot2::ggsave, ggsave_args)
}
write_table_list <- function(
    table_list,
    output_dir
) {
  for (table_name in names(table_list)) {
    table_data <- table_list[[table_name]]
    if (!is.data.frame(table_data)) {
      warning(
        "Output '",
        table_name,
        "' was not exported because it is not a table."
      )
      next
    }
    write.table(
      table_data,
      file.path(
        output_dir,
        paste0(table_name, ".tsv")
      ),
      sep = "\t",
      dec = ",",
      row.names = FALSE,
      quote = FALSE
    )
  }
}
tcr_list_to_table <- function(
    tcr_list,
    cell_table
) {
  required_cols <- c(
    "sample_id",
    "tcr"
  )
  missing_cols <- setdiff(
    required_cols,
    colnames(cell_table)
  )
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in 'cell_table': ",
      paste(missing_cols, collapse = ", ")
    )
  }
  tcr_table <- dplyr::bind_rows(
    lapply(
      names(tcr_list),
      function(condition) {
        tibble::tibble(
          condition = condition,
          tcr = as.character(
            tcr_list[[condition]]
          )
        )
      }
    )
  ) %>%
    dplyr::distinct()
  if (nrow(tcr_table) == 0L) {
    return(
      tibble::tibble(
        condition = character(),
        tcr = character(),
        tcr_cells = integer(),
        total_tcr_cells = integer(),
        tcr_frequency = numeric()
      )
    )
  }
  frequency_table <- cell_table %>%
    dplyr::mutate(
      sample_id = as.character(sample_id),
      tcr = as.character(tcr)
    ) %>%
    dplyr::count(
      sample_id,
      tcr,
      name = "tcr_cells"
    ) %>%
    dplyr::group_by(
      sample_id
    ) %>%
    dplyr::mutate(
      total_tcr_cells = sum(tcr_cells),
      tcr_frequency = 100 *
        tcr_cells /
        total_tcr_cells
    ) %>%
    dplyr::ungroup()
  tcr_table %>%
    dplyr::left_join(
      frequency_table,
      by = c(
        "condition" = "sample_id",
        "tcr" = "tcr"
      )
    ) %>%
    dplyr::arrange(
      condition,
      dplyr::desc(tcr_frequency),
      dplyr::desc(tcr_cells),
      tcr
    )
}
plot_TCR_clonotype_pies <- function(
    obj,
    analysis_conditions,
    compartment,
    coverage_targets = c(70, 100),
    max_clonotypes = 5000,
    legend_max = 20,
    save_freqs = NULL
) {
  required_cols <- c(
    "sample_id", "celltype", "tag.tcr_combi",
    "cdr3_combi", "v_gene_combi", "j_gene_combi",
    "chains_combi", "d_gene", "c_gene",
    "reads", "umis", "productive", "full_length", "paired"
  )
  missing_cols <- setdiff(
    required_cols,
    colnames(obj@meta.data)
  )
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in obj@meta.data: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  analysis_conditions <- unique(
    as.character(analysis_conditions)
  )
  coverage_targets <- unique(
    as.numeric(coverage_targets)
  )
  if (length(analysis_conditions) == 0) {
    stop(
      "'analysis_conditions' cannot be empty."
    )
  }
  if (
    length(coverage_targets) == 0 ||
    any(!is.finite(coverage_targets)) ||
    any(coverage_targets <= 0) ||
    any(coverage_targets > 100)
  ) {
    stop(
      "'coverage_targets' must contain values between 0 and 100."
    )
  }
  if (
    !is.numeric(max_clonotypes) ||
    length(max_clonotypes) != 1 ||
    max_clonotypes < 1
  ) {
    stop(
      "'max_clonotypes' must be a positive number."
    )
  }
  if (
    !is.numeric(legend_max) ||
    length(legend_max) != 1 ||
    legend_max < 1
  ) {
    stop(
      "'legend_max' must be a positive number."
    )
  }
  available_conditions <- obj@meta.data %>%
    dplyr::filter(
      celltype == compartment,
      tag.tcr_combi == "yes"
    ) %>%
    dplyr::pull(
      sample_id
    ) %>%
    unique() %>%
    as.character()
  missing_conditions <- setdiff(
    analysis_conditions,
    available_conditions
  )
  if (length(missing_conditions) > 0) {
    stop(
      "No TCR-positive cells were found for: ",
      paste(missing_conditions, collapse = ", ")
    )
  }
  if (!is.null(save_freqs)) {
    dir.create(
      save_freqs,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
  sample_names <- stats::setNames(
    analysis_conditions,
    analysis_conditions
  )
  freq_tables <- list()
  topx_lists <- list()
  plots <- list()
  topx_reports <- list()
  for (coverage_target in coverage_targets) {
    coverage_topx_lists <- list()
    for (i in names(sample_names)) {
      df <- obj@meta.data %>%
        dplyr::filter(
          celltype == compartment,
          tag.tcr_combi == "yes",
          sample_id == i
        ) %>%
        dplyr::mutate(
          cdr3_vj = paste(
            cdr3_combi,
            v_gene_combi,
            j_gene_combi,
            sep = "_"
          ),
          cdr3_vj = vapply(
            cdr3_vj,
            function(x) {
              if (length(x) == 0 || all(is.na(x))) {
                return("no_tcr")
              }
              if (is.list(x)) {
                x <- paste(
                  unlist(x),
                  collapse = ";"
                )
              }
              as.character(x)
            },
            character(1)
          )
        )
      df$cdr3_vj[
        is.na(df$cdr3_vj) |
          df$cdr3_vj == "NA_NA_NA"
      ] <- "no_tcr"
      sample_label <- sample_names[[i]]
      # --- Frequency table ---
      freq_table <- df %>%
        dplyr::filter(
          cdr3_vj != "no_tcr"
        ) %>%
        dplyr::group_by(
          sample_id,
          cdr3_vj
        ) %>%
        dplyr::summarise(
          n = dplyr::n(),
          chains_combi = dplyr::first(
            stats::na.omit(chains_combi),
            default = NA
          ),
          cdr3_combi = dplyr::first(
            stats::na.omit(cdr3_combi),
            default = NA
          ),
          v_gene_combi = dplyr::first(
            stats::na.omit(v_gene_combi),
            default = NA
          ),
          d_gene = dplyr::first(
            stats::na.omit(d_gene),
            default = NA
          ),
          j_gene_combi = dplyr::first(
            stats::na.omit(j_gene_combi),
            default = NA
          ),
          c_gene = dplyr::first(
            stats::na.omit(c_gene),
            default = NA
          ),
          reads = dplyr::first(
            stats::na.omit(reads),
            default = NA
          ),
          umis = dplyr::first(
            stats::na.omit(umis),
            default = NA
          ),
          productive = dplyr::first(
            stats::na.omit(productive),
            default = NA
          ),
          full_length = dplyr::first(
            stats::na.omit(full_length),
            default = NA
          ),
          paired = dplyr::first(
            stats::na.omit(paired),
            default = NA
          ),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          freq = n / sum(n) * 100
        ) %>%
        dplyr::arrange(
          dplyr::desc(freq)
        ) %>%
        dplyr::select(
          cdr3_vj,
          freq,
          n,
          chains_combi,
          cdr3_combi,
          v_gene_combi,
          d_gene,
          j_gene_combi,
          c_gene,
          reads,
          umis,
          productive,
          full_length,
          paired
        )
      if (nrow(freq_table) == 0) {
        warning(
          "No valid clonotypes were found for: ",
          sample_label
        )
        next
      }
      topx_tbl <- freq_table %>%
        dplyr::arrange(
          dplyr::desc(freq)
        ) %>%
        dplyr::mutate(
          cumfreq = cumsum(freq)
        )
      idx <- which(
        topx_tbl$cumfreq >= coverage_target
      )[1]
      if (is.na(idx)) {
        idx <- nrow(topx_tbl)
      }
      idx <- min(
        idx,
        max_clonotypes
      )
      topx_tbl <- topx_tbl %>%
        dplyr::slice(
          seq_len(idx)
        ) %>%
        dplyr::select(
          -cumfreq
        )
      full_table_name <- paste0(
        toupper(compartment),
        "_",
        sample_label,
        "_full"
      )
      top_table_name <- paste0(
        toupper(compartment),
        "_",
        sample_label,
        "_top",
        coverage_target
      )
      freq_tables[[full_table_name]] <- freq_table
      freq_tables[[top_table_name]] <- topx_tbl
      top_keys <- topx_tbl$cdr3_vj
      topx_table <- topx_tbl %>%
        dplyr::mutate(
          sample_label = sample_label
        ) %>%
        dplyr::select(
          sample_label,
          cdr3_vj,
          freq,
          n
        )
      topx_lists[[
        paste0(
          sample_label,
          "_top",
          coverage_target
        )
      ]] <- topx_table
      coverage_topx_lists[[
        sample_label
      ]] <- topx_table
      legend_keys <- head(
        top_keys,
        legend_max
      )
      inside_top_hidden_label <- "legend not shown"
      outside_top_label <- paste0(
        "outside top ",
        coverage_target,
        "%"
      )
      n_clones <- length(
        legend_keys
      )
      base_palette <- grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(
          8,
          "Paired"
        )
      )
      base_colors <- base_palette(
        max(
          11,
          n_clones
        )
      )
      color_manual <- stats::setNames(
        base_colors[
          seq_len(n_clones)
        ],
        legend_keys
      )
      color_manual[
        inside_top_hidden_label
      ] <- "gray70"
      color_manual[
        outside_top_label
      ] <- "gray90"
      cdr3_counts_sample <- freq_table %>%
        dplyr::mutate(
          cdr3_group = dplyr::case_when(
            cdr3_vj %in% legend_keys ~ cdr3_vj,
            cdr3_vj %in% top_keys ~ inside_top_hidden_label,
            TRUE ~ outside_top_label
          )
        ) %>%
        dplyr::group_by(
          cdr3_group
        ) %>%
        dplyr::summarise(
          n = sum(n),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          freq = n / sum(n) * 100
        ) %>%
        dplyr::arrange(
          dplyr::desc(freq)
        )
      plot_levels <- c(
        legend_keys[
          legend_keys %in%
            cdr3_counts_sample$cdr3_group
        ],
        inside_top_hidden_label,
        outside_top_label
      )
      plot_levels <- unique(
        plot_levels[
          plot_levels %in%
            cdr3_counts_sample$cdr3_group
        ]
      )
      cdr3_counts_sample <- cdr3_counts_sample %>%
        dplyr::mutate(
          cdr3_group = factor(
            cdr3_group,
            levels = plot_levels
          )
        )
      pie_sample <- ggplot2::ggplot(
        cdr3_counts_sample,
        ggplot2::aes(
          x = "",
          y = freq,
          fill = cdr3_group
        )
      ) +
        ggplot2::geom_col(
          color = "black",
          linewidth = 0.1
        ) +
        ggplot2::coord_polar(
          theta = "y",
          direction = 1
        ) +
        ggplot2::scale_fill_manual(
          values = color_manual,
          drop = FALSE
        ) +
        ggplot2::labs(
          title = sample_label,
          fill = paste0(
            "CDR3_VJ (top ",
            coverage_target,
            "%; max ",
            max_clonotypes,
            ")"
          )
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(
          legend.position = "right",
          plot.title = ggplot2::element_text(
            hjust = 0,
            margin = ggplot2::margin(
              l = 40
            ),
            size = 10
          ),
          plot.margin = ggplot2::margin(
            t = 10,
            r = 5,
            b = 5,
            l = 5
          )
        )
      plot_name <- paste0(
        toupper(compartment),
        "_",
        sample_label,
        "_top",
        coverage_target
      )
      plots[[plot_name]] <- pie_sample
    }
    topx_report <- dplyr::bind_rows(
      coverage_topx_lists
    ) %>%
      dplyr::select(
        sample_label,
        cdr3_vj
      ) %>%
      dplyr::group_by(
        cdr3_vj
      ) %>%
      dplyr::summarise(
        sample_label = paste(
          sort(
            unique(
              stats::na.omit(
                sample_label
              )
            )
          ),
          collapse = " | "
        ),
        .groups = "drop"
      )
    report_name <- paste0(
      toupper(compartment),
      "_TCR_interrogated_top",
      coverage_target
    )
    topx_reports[[report_name]] <- topx_report
    if (!is.null(save_freqs)) {
      write.table(
        topx_report,
        file.path(
          save_freqs,
          paste0(
            report_name,
            ".tsv"
          )
        ),
        sep = "\t",
        dec = ",",
        row.names = FALSE,
        quote = FALSE
      )
    }
  }
  return(
    list(
      plots = plots,
      freq_tables = freq_tables,
      topx_lists = topx_lists,
      topx_reports = topx_reports,
      compartment = compartment,
      analysis_conditions = analysis_conditions,
      coverage_targets = coverage_targets,
      max_clonotypes = max_clonotypes,
      legend_max = legend_max
    )
  )
}
plot_global_TCR_intersections <- function(
    obj,
    compartment,
    conditions = NULL
) {
  required_cols <- c(
    "sample_id", "celltype", "tag.tcr_combi",
    "cdr3_combi", "v_gene_combi", "j_gene_combi"
  )
  missing_cols <- setdiff(
    required_cols,
    colnames(obj@meta.data)
  )
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in obj@meta.data: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
    stop(
      "Package 'ggVennDiagram' is required. Install it with: ",
      "install.packages('ggVennDiagram')"
    )
  }
  df <- obj@meta.data %>%
    tibble::rownames_to_column(
      var = "cell_id"
    ) %>%
    dplyr::filter(
      celltype == compartment,
      tag.tcr_combi == "yes",
      !is.na(sample_id),
      !is.na(cdr3_combi),
      !is.na(v_gene_combi),
      !is.na(j_gene_combi)
    ) %>%
    dplyr::mutate(
      sample_id = as.character(sample_id),
      tcr = paste(
        cdr3_combi,
        v_gene_combi,
        j_gene_combi,
        sep = "_"
      )
    ) %>%
    dplyr::filter(
      tcr != "NA_NA_NA"
    )
  if (nrow(df) == 0) {
    stop(
      "No valid TCR-positive cells were found for compartment: ",
      compartment
    )
  }
  available_conditions <- unique(
    df$sample_id
  )
  if (is.null(conditions)) {
    conditions <- available_conditions
  }
  conditions <- unique(
    as.character(conditions)
  )
  missing_conditions <- setdiff(
    conditions,
    available_conditions
  )
  if (length(missing_conditions) > 0) {
    stop(
      "No valid TCR-positive cells were found for: ",
      paste(missing_conditions, collapse = ", ")
    )
  }
  if (length(conditions) < 2) {
    stop(
      "At least two conditions are required."
    )
  }
  if (length(conditions) > 4) {
    stop(
      "The global Venn diagram supports a maximum of four conditions. ",
      "Use an UpSet plot for more than four conditions."
    )
  }
  df <- df %>%
    dplyr::filter(
      sample_id %in% conditions
    )
  tcr_presence_long <- df %>%
    dplyr::select(
      sample_id,
      tcr
    ) %>%
    dplyr::distinct()
  tcr_sets_table <- tcr_presence_long %>%
    dplyr::group_by(
      sample_id
    ) %>%
    dplyr::summarise(
      tcrs = list(
        unique(tcr)
      ),
      .groups = "drop"
    )
  tcr_sets <- stats::setNames(
    tcr_sets_table$tcrs,
    tcr_sets_table$sample_id
  )
  tcr_sets <- tcr_sets[
    conditions
  ]
  membership_table <- tcr_presence_long %>%
    dplyr::mutate(
      present = TRUE
    ) %>%
    tidyr::pivot_wider(
      names_from = sample_id,
      values_from = present,
      values_fill = FALSE
    ) %>%
    dplyr::select(
      tcr,
      dplyr::all_of(conditions)
    )
  membership_matrix <- as.matrix(
    membership_table[
      conditions
    ]
  )
  storage.mode(membership_matrix) <- "logical"
  intersection <- apply(
    membership_matrix,
    1,
    function(present) {
      selected_conditions <- conditions[
        present
      ]
      if (length(selected_conditions) == 1) {
        paste0(
          selected_conditions,
          " only"
        )
      } else {
        paste(
          selected_conditions,
          collapse = " + "
        )
      }
    }
  )
  intersection_key <- apply(
    membership_matrix,
    1,
    function(present) {
      paste(
        conditions[present],
        collapse = "__"
      )
    }
  )
  membership_table <- membership_table %>%
    dplyr::mutate(
      n_conditions = rowSums(
        membership_matrix
      ),
      intersection = intersection,
      intersection_key = intersection_key
    )
  clonotype_table <- membership_table %>%
    dplyr::select(
      tcr,
      intersection,
      intersection_key,
      n_conditions,
      dplyr::all_of(conditions)
    ) %>%
    dplyr::arrange(
      dplyr::desc(n_conditions),
      intersection,
      tcr
    )
  summary_table <- clonotype_table %>%
    dplyr::count(
      intersection,
      intersection_key,
      n_conditions,
      name = "n_clonotypes"
    ) %>%
    dplyr::arrange(
      dplyr::desc(n_conditions),
      dplyr::desc(n_clonotypes),
      intersection
    )
  intersection_tcrs <- split(
    clonotype_table$tcr,
    clonotype_table$intersection
  )
  cell_table <- df %>%
    dplyr::left_join(
      clonotype_table %>%
        dplyr::select(
          tcr,
          intersection,
          intersection_key,
          n_conditions
        ),
      by = "tcr"
    )
  cell_summary_table <- cell_table %>%
    dplyr::count(
      intersection,
      intersection_key,
      sample_id,
      name = "n_cells"
    ) %>%
    dplyr::mutate(
      sample_id = factor(
        sample_id,
        levels = conditions
      )
    ) %>%
    dplyr::arrange(
      intersection,
      sample_id
    )
  venn_plot <- ggVennDiagram::ggVennDiagram(
    x = tcr_sets,
    label = "count",
    label_alpha = 0,
    edge_size = 0.6,
    set_size = 4,
    label_size = 3.5
  ) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "grey70"
    ) +
    ggplot2::coord_equal(
      clip = "off"
    ) +
    ggplot2::labs(
      title = paste0(
        toupper(compartment),
        ": global TCR overlap"
      ),
      subtitle = paste0(
        nrow(clonotype_table),
        " total clonotypes"
      )
    ) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 14
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        size = 10
      ),
      plot.margin = ggplot2::margin(
        t = 30,
        r = 100,
        b = 30,
        l = 100
      )
    )
  plot_name <- paste0(
    toupper(compartment),
    "_global_TCR_Venn"
  )
  return(
    list(
      plot = venn_plot,
      plot_name = plot_name,
      compartment = compartment,
      conditions = conditions,
      tcr_sets = tcr_sets,
      membership_table = membership_table,
      clonotype_table = clonotype_table,
      summary_table = summary_table,
      intersection_tcrs = intersection_tcrs,
      cell_table = cell_table,
      cell_summary_table = cell_summary_table
    )
  )
}
plot_TCR_only_vs_non_only_marker_expression <- function(
    global_tcr_result,
    selected_conditions,
    markers,
    control_condition = "Unpulsed",
    show_mean = TRUE,
    ncol = 6
) {
  required_elements <- c(
    "compartment", "conditions", "clonotype_table",
    "cell_table", "intersection_tcrs"
  )
  missing_elements <- setdiff(
    required_elements,
    names(global_tcr_result)
  )
  if (length(missing_elements) > 0) {
    stop(
      "Missing elements in 'global_tcr_result': ",
      paste(missing_elements, collapse = ", ")
    )
  }
  selected_conditions <- unique(
    as.character(selected_conditions)
  )
  markers <- unique(
    as.character(markers)
  )
  if (length(selected_conditions) == 0) {
    stop(
      "'selected_conditions' cannot be empty."
    )
  }
  if (length(markers) == 0) {
    stop(
      "'markers' cannot be empty."
    )
  }
  if (length(control_condition) != 1) {
    stop(
      "'control_condition' must contain exactly one condition."
    )
  }
  invalid_conditions <- setdiff(
    selected_conditions,
    global_tcr_result$conditions
  )
  if (length(invalid_conditions) > 0) {
    stop(
      "Conditions not found in the global analysis: ",
      paste(invalid_conditions, collapse = ", ")
    )
  }
  if (!control_condition %in% global_tcr_result$conditions) {
    stop(
      "Control condition not found in the global analysis: ",
      control_condition
    )
  }
  if (control_condition %in% selected_conditions) {
    stop(
      "'control_condition' cannot also be included in 'selected_conditions'."
    )
  }
  marker_cols <- paste0(
    markers,
    "_mag"
  )
  missing_marker_cols <- setdiff(
    marker_cols,
    colnames(global_tcr_result$cell_table)
  )
  if (length(missing_marker_cols) > 0) {
    stop(
      "Missing marker columns in 'global_tcr_result$cell_table': ",
      paste(missing_marker_cols, collapse = ", ")
    )
  }
  control_group <- paste0(
    control_condition,
    " all TCRs"
  )
  control_data <- global_tcr_result$cell_table %>%
    dplyr::filter(
      sample_id == control_condition
    ) %>%
    dplyr::mutate(
      comparison_class = "control",
      plot_group = control_group
    )
  if (nrow(control_data) == 0) {
    stop(
      "No TCR-positive cells were found for the control condition: ",
      control_condition
    )
  }
  condition_data_list <- list()
  only_tcrs <- list()
  non_only_tcrs <- list()
  for (condition in selected_conditions) {
    only_intersection <- paste0(
      condition,
      " only"
    )
    condition_data <- global_tcr_result$cell_table %>%
      dplyr::filter(
        sample_id == condition
      )
    if (nrow(condition_data) == 0) {
      stop(
        "No TCR-positive cells were found for condition: ",
        condition
      )
    }
    condition_only_data <- condition_data %>%
      dplyr::filter(
        intersection == only_intersection
      ) %>%
      dplyr::mutate(
        comparison_class = "only",
        plot_group = paste0(
          condition,
          " only"
        )
      )
    condition_non_only_data <- condition_data %>%
      dplyr::filter(
        intersection != only_intersection
      ) %>%
      dplyr::mutate(
        comparison_class = "non-only",
        plot_group = paste0(
          condition,
          " non-only"
        )
      )
    if (nrow(condition_only_data) == 0) {
      warning(
        "No cells were found for: ",
        only_intersection
      )
    }
    if (nrow(condition_non_only_data) == 0) {
      warning(
        "No non-only cells were found for condition: ",
        condition
      )
    }
    only_tcrs[[condition]] <- unique(
      condition_only_data$tcr
    )
    non_only_tcrs[[condition]] <- unique(
      condition_non_only_data$tcr
    )
    condition_data_list[[condition]] <- dplyr::bind_rows(
      condition_only_data,
      condition_non_only_data
    )
  }
  queried_data <- dplyr::bind_rows(
    condition_data_list
  )
  plot_data <- dplyr::bind_rows(
    control_data,
    queried_data
  )
  group_levels <- control_group
  for (condition in selected_conditions) {
    group_levels <- c(
      group_levels,
      paste0(condition, " only"),
      paste0(condition, " non-only")
    )
  }
  plot_data <- plot_data %>%
    dplyr::mutate(
      plot_group = factor(
        plot_group,
        levels = group_levels
      ),
      sample_id = factor(
        sample_id,
        levels = c(
          control_condition,
          selected_conditions
        )
      ),
      comparison_class = factor(
        comparison_class,
        levels = c(
          "control",
          "only",
          "non-only"
        )
      )
    )
  selected_condition_colors <- stats::setNames(
    grDevices::hcl.colors(
      length(selected_conditions),
      palette = "Dark 3"
    ),
    selected_conditions
  )
  condition_colors <- c(
    stats::setNames(
      "grey70",
      control_condition
    ),
    selected_condition_colors
  )
  group_summary <- plot_data %>%
    dplyr::group_by(
      plot_group,
      sample_id,
      comparison_class
    ) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      n_clonotypes = dplyr::n_distinct(
        tcr
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      plot_group
    )
  non_only_intersection_summary <- queried_data %>%
    dplyr::filter(
      comparison_class == "non-only"
    ) %>%
    dplyr::group_by(
      sample_id,
      intersection
    ) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      n_clonotypes = dplyr::n_distinct(
        tcr
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      sample_id,
      dplyr::desc(n_clonotypes),
      intersection
    )
  marker_plots <- list()
  marker_summary_tables <- list()
  for (marker in markers) {
    marker_col <- paste0(
      marker,
      "_mag"
    )
    marker_plot_data <- plot_data %>%
      dplyr::transmute(
        sample_id = sample_id,
        intersection = intersection,
        comparison_class = comparison_class,
        plot_group = plot_group,
        tcr = tcr,
        expression = .data[[marker_col]]
      )
    marker_summary_tables[[marker]] <- marker_plot_data %>%
      dplyr::group_by(
        plot_group,
        sample_id,
        comparison_class
      ) %>%
      dplyr::summarise(
        marker = marker,
        n_cells = sum(
          !is.na(expression)
        ),
        n_clonotypes = dplyr::n_distinct(
          tcr
        ),
        mean_expression = mean(
          expression,
          na.rm = TRUE
        ),
        median_expression = stats::median(
          expression,
          na.rm = TRUE
        ),
        .groups = "drop"
      )
    marker_plot <- ggplot2::ggplot(
      marker_plot_data,
      ggplot2::aes(
        x = plot_group,
        y = expression,
        fill = sample_id
      )
    ) +
      ggplot2::geom_boxplot(
        outlier.shape = NA,
        linewidth = 0.3,
        width = 0.7
      )
    if (show_mean) {
      marker_plot <- marker_plot +
        ggplot2::stat_summary(
          fun = mean,
          geom = "point",
          shape = 23,
          size = 2,
          fill = "white",
          na.rm = TRUE
        )
    }
    marker_plots[[marker]] <- marker_plot +
      ggplot2::scale_x_discrete(
        limits = group_levels,
        drop = FALSE
      ) +
      ggplot2::scale_fill_manual(
        values = condition_colors,
        drop = FALSE
      ) +
      ggplot2::labs(
        title = marker,
        x = NULL,
        y = paste0(
          marker,
          " expression"
        ),
        fill = "Condition"
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(
          hjust = 0.5
        ),
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1,
          vjust = 1
        ),
        axis.line = ggplot2::element_line(
          linewidth = 0.3
        ),
        axis.ticks = ggplot2::element_line(
          linewidth = 0.3
        ),
        aspect.ratio = 1
      )
  }
  final_plot <- patchwork::wrap_plots(
    marker_plots[markers],
    ncol = ncol,
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = paste0(
        toupper(global_tcr_result$compartment),
        ": condition-specific vs shared clonotypes"
      ),
      subtitle = paste0(
        paste(
          selected_conditions,
          collapse = " | "
        ),
        " vs ",
        control_condition,
        " all TCRs"
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold"
        ),
        plot.subtitle = ggplot2::element_text(
          hjust = 0.5
        )
      )
    ) &
    ggplot2::theme(
      legend.position = "bottom"
    )
  marker_summary_table <- dplyr::bind_rows(
    marker_summary_tables
  )
  plot_name <- paste0(
    toupper(global_tcr_result$compartment),
    "_",
    paste(
      selected_conditions,
      collapse = "_"
    ),
    "_only_vs_non_only"
  )
  return(
    list(
      plot = final_plot,
      plot_name = plot_name,
      selected_conditions = selected_conditions,
      control_condition = control_condition,
      plot_data = plot_data,
      control_data = control_data,
      queried_data = queried_data,
      only_tcrs = only_tcrs,
      non_only_tcrs = non_only_tcrs,
      group_summary = group_summary,
      non_only_intersection_summary = non_only_intersection_summary,
      marker_summary_table = marker_summary_table,
      marker_plots = marker_plots
    )
  )
}
plot_TCR_frequency_FC_marker_expression <- function(
    obj,
    selected_conditions,
    markers,
    control_condition = "Unpulsed",
    compartment,
    frequency_threshold = 5,
    fc_threshold = 2,
    ncol = 6
) {
  selected_conditions <- unique(
    as.character(selected_conditions)
  )
  markers <- unique(
    as.character(markers)
  )
  if (length(selected_conditions) == 0) {
    stop(
      "'selected_conditions' cannot be empty."
    )
  }
  if (length(markers) == 0) {
    stop(
      "'markers' cannot be empty."
    )
  }
  if (length(control_condition) != 1) {
    stop(
      "'control_condition' must contain exactly one condition."
    )
  }
  if (control_condition %in% selected_conditions) {
    stop(
      "'control_condition' cannot also be included in 'selected_conditions'."
    )
  }
  if (
    !is.numeric(frequency_threshold) ||
    length(frequency_threshold) != 1 ||
    frequency_threshold < 0
  ) {
    stop(
      "'frequency_threshold' must be a single non-negative number."
    )
  }
  if (
    !is.numeric(fc_threshold) ||
    length(fc_threshold) != 1 ||
    fc_threshold <= 0
  ) {
    stop(
      "'fc_threshold' must be a single positive number."
    )
  }
  marker_cols <- paste0(
    markers,
    "_mag"
  )
  required_cols <- c(
    "sample_id", "celltype", "tag.tcr_combi",
    "cdr3_combi", "v_gene_combi", "j_gene_combi",
    marker_cols
  )
  missing_cols <- setdiff(
    required_cols,
    colnames(obj@meta.data)
  )
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in obj@meta.data: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  conditions <- c(
    control_condition,
    selected_conditions
  )
  df <- obj@meta.data %>%
    dplyr::filter(
      sample_id %in% conditions,
      celltype == compartment,
      tag.tcr_combi == "yes",
      !is.na(sample_id),
      !is.na(cdr3_combi),
      !is.na(v_gene_combi),
      !is.na(j_gene_combi)
    ) %>%
    dplyr::mutate(
      sample_id = as.character(sample_id),
      tcr = paste(
        cdr3_combi,
        v_gene_combi,
        j_gene_combi,
        sep = "_"
      )
    ) %>%
    dplyr::filter(
      tcr != "NA_NA_NA"
    )
  if (nrow(df) == 0) {
    stop(
      "No valid TCR-positive cells were found for the selected compartment ",
      "and conditions."
    )
  }
  available_conditions <- unique(
    df$sample_id
  )
  missing_conditions <- setdiff(
    conditions,
    available_conditions
  )
  if (length(missing_conditions) > 0) {
    stop(
      "No valid TCR-positive cells were found for: ",
      paste(missing_conditions, collapse = ", ")
    )
  }
  group_totals <- df %>%
    dplyr::count(
      sample_id,
      name = "total_tcr_cells"
    )
  frequency_table <- df %>%
    dplyr::count(
      sample_id,
      tcr,
      name = "tcr_cells"
    ) %>%
    dplyr::left_join(
      group_totals,
      by = "sample_id"
    ) %>%
    dplyr::mutate(
      tcr_frequency = 100 *
        tcr_cells /
        total_tcr_cells
    ) %>%
    dplyr::select(
      sample_id,
      tcr,
      tcr_cells,
      total_tcr_cells,
      tcr_frequency
    )
  condition_colors <- stats::setNames(
    c(
      "grey70",
      grDevices::hcl.colors(
        length(selected_conditions),
        palette = "Dark 3"
      )
    ),
    conditions
  )
  comparison_results <- list()
  comparison_plots <- list()
  prioritization_tables <- list()
  frequency_tables <- list()
  for (condition in selected_conditions) {
    comparison_name <- paste0(
      control_condition,
      "_vs_",
      condition
    )
    control_tcrs <- unique(
      df$tcr[
        df$sample_id == control_condition
      ]
    )
    condition_tcrs <- unique(
      df$tcr[
        df$sample_id == condition
      ]
    )
    shared_tcrs <- intersect(
      control_tcrs,
      condition_tcrs
    )
    if (length(shared_tcrs) == 0) {
      warning(
        "No shared TCRs were found between ",
        control_condition,
        " and ",
        condition,
        "."
      )
      next
    }
    control_frequency_col <- paste0(
      "tcr_frequency_",
      control_condition
    )
    condition_frequency_col <- paste0(
      "tcr_frequency_",
      condition
    )
    control_cells_col <- paste0(
      "tcr_cells_",
      control_condition
    )
    condition_cells_col <- paste0(
      "tcr_cells_",
      condition
    )
    control_mean_col <- paste0(
      "mean_expression_",
      control_condition
    )
    condition_mean_col <- paste0(
      "mean_expression_",
      condition
    )
    comparison_frequency_table <- frequency_table %>%
      dplyr::filter(
        sample_id %in% c(
          control_condition,
          condition
        ),
        tcr %in% shared_tcrs
      ) %>%
      dplyr::select(
        sample_id,
        tcr,
        tcr_cells,
        tcr_frequency
      ) %>%
      tidyr::pivot_wider(
        names_from = sample_id,
        values_from = c(
          tcr_cells,
          tcr_frequency
        )
      ) %>%
      dplyr::mutate(
        frequency_increment =
          .data[[condition_frequency_col]] -
          .data[[control_frequency_col]]
      ) %>%
      dplyr::filter(
        !is.na(frequency_increment),
        is.finite(frequency_increment)
      ) %>%
      dplyr::arrange(
        dplyr::desc(frequency_increment),
        dplyr::desc(
          .data[[condition_frequency_col]]
        )
      )
    frequency_tables[[comparison_name]] <- comparison_frequency_table
    comparison_marker_tables <- list()
    comparison_marker_plots <- list()
    for (marker in markers) {
      marker_col <- paste0(
        marker,
        "_mag"
      )
      marker_positive_values <- df[[marker_col]][
        df$sample_id %in% c(
          control_condition,
          condition
        ) &
          df$tcr %in% shared_tcrs &
          is.finite(df[[marker_col]]) &
          df[[marker_col]] > 0
      ]
      if (length(marker_positive_values) == 0) {
        warning(
          "No positive expression values were found for marker ",
          marker,
          " in ",
          comparison_name,
          "."
        )
        next
      }
      marker_min_positive <- min(
        marker_positive_values,
        na.rm = TRUE
      )
      expression_table <- df %>%
        dplyr::filter(
          sample_id %in% c(
            control_condition,
            condition
          ),
          tcr %in% shared_tcrs
        ) %>%
        dplyr::group_by(
          tcr,
          sample_id
        ) %>%
        dplyr::summarise(
          mean_expression = mean(
            .data[[marker_col]],
            na.rm = TRUE
          ),
          .groups = "drop"
        ) %>%
        tidyr::pivot_wider(
          names_from = sample_id,
          values_from = mean_expression,
          names_prefix = "mean_expression_"
        ) %>%
        dplyr::mutate(
          "{control_mean_col}" := dplyr::case_when(
            is.na(.data[[control_mean_col]]) ~ NA_real_,
            !is.finite(.data[[control_mean_col]]) ~ NA_real_,
            .data[[control_mean_col]] <= 0 ~ marker_min_positive,
            TRUE ~ .data[[control_mean_col]]
          ),
          "{condition_mean_col}" := dplyr::case_when(
            is.na(.data[[condition_mean_col]]) ~ NA_real_,
            !is.finite(.data[[condition_mean_col]]) ~ NA_real_,
            .data[[condition_mean_col]] <= 0 ~ marker_min_positive,
            TRUE ~ .data[[condition_mean_col]]
          ),
          FC =
            .data[[condition_mean_col]] /
            .data[[control_mean_col]]
        )
      prioritization_table <- expression_table %>%
        dplyr::left_join(
          comparison_frequency_table,
          by = "tcr"
        ) %>%
        dplyr::filter(
          !is.na(FC),
          is.finite(FC),
          !is.na(frequency_increment),
          frequency_increment >= frequency_threshold,
          FC > fc_threshold
        ) %>%
        dplyr::transmute(
          condition = condition,
          control_condition = control_condition,
          marker = marker,
          tcr = tcr,
          mean_expression_control =
            .data[[control_mean_col]],
          mean_expression_condition =
            .data[[condition_mean_col]],
          tcr_cells_control =
            .data[[control_cells_col]],
          tcr_cells_condition =
            .data[[condition_cells_col]],
          tcr_frequency_control =
            .data[[control_frequency_col]],
          tcr_frequency_condition =
            .data[[condition_frequency_col]],
          frequency_increment =
            frequency_increment,
          FC = FC,
          frequency_threshold =
            frequency_threshold,
          fc_threshold =
            fc_threshold
        ) %>%
        dplyr::arrange(
          dplyr::desc(frequency_increment),
          dplyr::desc(FC),
          dplyr::desc(tcr_frequency_condition)
        )
      comparison_marker_tables[[marker]] <- prioritization_table
      selected_tcrs <- unique(
        prioritization_table$tcr
      )
      if (length(selected_tcrs) == 0L) {
        next
      }
      marker_plot_data <- df %>%
        dplyr::filter(
          sample_id %in% c(
            control_condition,
            condition
          ),
          tcr %in% selected_tcrs
        ) %>%
        dplyr::transmute(
          sample_id = factor(
            sample_id,
            levels = c(
              control_condition,
              condition
            )
          ),
          tcr = tcr,
          expression = .data[[marker_col]]
        )
      comparison_marker_plots[[marker]] <- ggplot2::ggplot(
        marker_plot_data,
        ggplot2::aes(
          x = sample_id,
          y = expression,
          fill = sample_id
        )
      ) +
        ggplot2::geom_boxplot(
          outlier.shape = NA,
          linewidth = 0.3,
          width = 0.7
        ) +
        ggplot2::stat_summary(
          fun = mean,
          geom = "point",
          shape = 23,
          size = 2,
          fill = "white",
          na.rm = TRUE
        ) +
        ggplot2::scale_x_discrete(
          limits = c(
            control_condition,
            condition
          ),
          drop = FALSE
        ) +
        ggplot2::scale_fill_manual(
          values = condition_colors[
            c(
              control_condition,
              condition
            )
          ],
          drop = FALSE
        ) +
        ggplot2::labs(
          title = marker,
          subtitle = paste0(
            length(selected_tcrs),
            " prioritized TCRs"
          ),
          x = NULL,
          y = paste0(
            marker,
            " expression"
          )
        ) +
        ggplot2::theme_classic() +
        ggplot2::theme(
          legend.position = "none",
          plot.title = ggplot2::element_text(
            hjust = 0.5
          ),
          plot.subtitle = ggplot2::element_text(
            hjust = 0.5
          ),
          axis.text.x = ggplot2::element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
          ),
          axis.line = ggplot2::element_line(
            linewidth = 0.3
          ),
          axis.ticks = ggplot2::element_line(
            linewidth = 0.3
          ),
          aspect.ratio = 1
        )
    }
    comparison_prioritization_table <- dplyr::bind_rows(
      comparison_marker_tables
    )
    prioritization_tables[[comparison_name]] <-
      comparison_prioritization_table
    available_markers <- markers[
      markers %in% names(
        comparison_marker_plots
      )
    ]
    if (length(available_markers) > 0) {
      comparison_plot <- patchwork::wrap_plots(
        comparison_marker_plots[
          available_markers
        ],
        ncol = ncol
      ) +
        patchwork::plot_annotation(
          title = paste0(
            toupper(compartment),
            ": ",
            condition,
            " vs ",
            control_condition
          ),
          subtitle = paste0(
            "Delta frequency >= ",
            frequency_threshold,
            " pp | expression FC > ",
            fc_threshold
          ),
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(
              hjust = 0.5,
              face = "bold"
            ),
            plot.subtitle = ggplot2::element_text(
              hjust = 0.5
            )
          )
        )
    } else {
      comparison_plot <- NULL
    }
    comparison_results[[comparison_name]] <- list(
      condition = condition,
      control_condition = control_condition,
      shared_tcrs = shared_tcrs,
      frequency_table = comparison_frequency_table,
      prioritization_table =
        comparison_prioritization_table,
      marker_tables = comparison_marker_tables,
      marker_plots = comparison_marker_plots,
      plot = comparison_plot
    )
    if (!is.null(comparison_plot)) {
      comparison_plots[[comparison_name]] <- comparison_plot
    }
  }
  prioritization_table <- dplyr::bind_rows(
    prioritization_tables
  )
  if (nrow(prioritization_table) > 0) {
    prioritization_table <- prioritization_table %>%
      dplyr::distinct(
        condition,
        marker,
        tcr,
        .keep_all = TRUE
      ) %>%
      dplyr::arrange(
        condition,
        dplyr::desc(frequency_increment),
        dplyr::desc(FC),
        marker,
        tcr
      )
  }
  if (nrow(prioritization_table) > 0) {
    prioritized_tcr_summary <- prioritization_table %>%
      dplyr::group_by(
        condition,
        control_condition,
        tcr
      ) %>%
      dplyr::summarise(
        responsible_markers = paste(
          sort(
            unique(marker)
          ),
          collapse = "; "
        ),
        n_responsible_markers =
          dplyr::n_distinct(marker),
        max_FC = max(
          FC,
          na.rm = TRUE
        ),
        mean_FC = mean(
          FC,
          na.rm = TRUE
        ),
        tcr_frequency_control =
          dplyr::first(
            tcr_frequency_control
          ),
        tcr_frequency_condition =
          dplyr::first(
            tcr_frequency_condition
          ),
        frequency_increment =
          dplyr::first(
            frequency_increment
          ),
        .groups = "drop"
      ) %>%
      dplyr::arrange(
        condition,
        dplyr::desc(frequency_increment),
        dplyr::desc(n_responsible_markers),
        dplyr::desc(max_FC)
      )
  } else {
    prioritized_tcr_summary <- tibble::tibble(
      condition = character(),
      control_condition = character(),
      tcr = character(),
      responsible_markers = character(),
      n_responsible_markers = integer(),
      max_FC = numeric(),
      mean_FC = numeric(),
      tcr_frequency_control = numeric(),
      tcr_frequency_condition = numeric(),
      frequency_increment = numeric()
    )
  }
  if (length(comparison_plots) > 0) {
    final_plot <- patchwork::wrap_plots(
      comparison_plots,
      ncol = 1
    )
  } else {
    final_plot <- NULL
  }
  plot_name <- paste0(
    toupper(compartment),
    "_",
    paste(
      selected_conditions,
      collapse = "_"
    ),
    "_frequency_FC_prioritization"
  )
  return(
    list(
      plot = final_plot,
      plot_name = plot_name,
      control_condition = control_condition,
      selected_conditions = selected_conditions,
      compartment = compartment,
      frequency_threshold = frequency_threshold,
      fc_threshold = fc_threshold,
      comparison_results = comparison_results,
      frequency_tables = frequency_tables,
      prioritization_table = prioritization_table,
      prioritized_tcr_summary = prioritized_tcr_summary
    )
  )
}
generate_tcr_rank_plots <- function(
    obj_sub,
    sample_id_value,
    genes,
    compartment,
    region_size = 20L,
    gene_colors = NULL,
    minimum_markers_to_color = 2L,
    minimum_markers_to_label = 3L,
    annotation_y = 0.100,
    gene_width = 4,
    barplot_ncol = 2L,
    umap_background_size = 0.25,
    umap_selected_size = 0.75
) {
  genes_expression <- paste(
    deparse(substitute(genes)),
    collapse = ""
  )
  calling_environment <- parent.frame()
  genes <- unique(as.character(genes))
  is_positive_set <- exists(
    "genes_magic_pos",
    envir = calling_environment,
    inherits = TRUE
  ) && setequal(
    genes,
    get(
      "genes_magic_pos",
      envir = calling_environment,
      inherits = TRUE
    )
  )
  is_negative_set <- exists(
    "genes_magic_neg",
    envir = calling_environment,
    inherits = TRUE
  ) && setequal(
    genes,
    get(
      "genes_magic_neg",
      envir = calling_environment,
      inherits = TRUE
    )
  )
  gene_set_name <- if (
    genes_expression == "genes_magic_pos" ||
    is_positive_set
  ) {
    "genes_magic_pos"
  } else if (
    genes_expression == "genes_magic_neg" ||
    is_negative_set
  ) {
    "genes_magic_neg"
  } else {
    "custom_genes"
  }
  if (is.null(gene_colors)) {
    color_table_name <- switch(
      gene_set_name,
      genes_magic_pos = "top_genes",
      genes_magic_neg = "negative_control_genes",
      NULL
    )
    if (
      !is.null(color_table_name) &&
      exists(
        color_table_name,
        envir = calling_environment,
        inherits = TRUE
      )
    ) {
      color_table <- get(
        color_table_name,
        envir = calling_environment,
        inherits = TRUE
      )
      gene_colors <- setNames(
        color_table$color,
        color_table$gene
      )
    } else {
      gene_colors <- setNames(
        pals::glasbey(length(genes)),
        genes
      )
    }
  }
  missing_gene_colors <- setdiff(
    genes,
    names(gene_colors)
  )
  if (length(missing_gene_colors)) {
    stop(
      "Missing colors for: ",
      paste(missing_gene_colors, collapse = ", ")
    )
  }
  gene_colors <- gene_colors[genes]
  region_size <- as.integer(region_size)
  minimum_markers_to_color <- as.integer(
    minimum_markers_to_color
  )
  minimum_markers_to_label <- as.integer(
    minimum_markers_to_label
  )
  if (!length(genes)) {
    stop("The genes argument is empty.")
  }
  if (
    length(region_size) != 1L ||
    is.na(region_size) ||
    region_size < 1L
  ) {
    stop("region_size must be a positive integer.")
  }
  if (
    !is.numeric(umap_background_size) ||
    !is.numeric(umap_selected_size) ||
    umap_background_size <= 0 ||
    umap_selected_size <= umap_background_size
  ) {
    stop(
      "umap_selected_size must be greater than umap_background_size."
    )
  }
  required_metadata <- c(
    "sample_id",
    "TCR_ID",
    "CloneID"
  )
  missing_metadata <- setdiff(
    required_metadata,
    colnames(obj_sub[[]])
  )
  if (length(missing_metadata)) {
    stop(
      "Missing metadata columns: ",
      paste(missing_metadata, collapse = ", ")
    )
  }
  if (!sample_id_value %in% obj_sub$sample_id) {
    stop(
      "sample_id '",
      sample_id_value,
      "' was not found."
    )
  }
  magic_columns <- paste0(
    genes,
    "_mag"
  )
  missing_magic_columns <- setdiff(
    magic_columns,
    colnames(obj_sub[[]])
  )
  if (length(missing_magic_columns)) {
    stop(
      "Missing MAGIC metadata columns: ",
      paste(missing_magic_columns, collapse = ", ")
    )
  }
  sample_cells <- rownames(obj_sub[[]])[
    !is.na(obj_sub$sample_id) &
      obj_sub$sample_id == sample_id_value
  ]
  obj_sample <- subset(
    obj_sub,
    cells = sample_cells
  )
  valid_clone <- !is.na(obj_sample$TCR_ID) &
    !is.na(obj_sample$CloneID)
  if (!any(valid_clone)) {
    stop(
      "No valid clonotypes were found for ",
      sample_id_value,
      "."
    )
  }
  clone_freq <- table(
    obj_sample$CloneID[valid_clone]
  )
  obj_sample$CloneFrequency <- NA_integer_
  obj_sample$CloneFrequency[valid_clone] <- as.integer(
    clone_freq[obj_sample$CloneID[valid_clone]]
  )
  obj_sample$CloneProportion <- NA_real_
  obj_sample$CloneProportion[valid_clone] <-
    obj_sample$CloneFrequency[valid_clone] /
    sum(valid_clone)
  clone_reference <- obj_sample[[]] |>
    dplyr::filter(
      !is.na(TCR_ID),
      !is.na(CloneID)
    ) |>
    dplyr::group_by(
      TCR_ID,
      CloneID
    ) |>
    dplyr::summarise(
      CloneFrequency = dplyr::first(CloneFrequency),
      CloneProportion = dplyr::first(CloneProportion),
      .groups = "drop"
    )
  clone_proportion_sum <- sum(
    clone_reference$CloneProportion,
    na.rm = TRUE
  )
  stopifnot(
    isTRUE(all.equal(
      clone_proportion_sum,
      1,
      tolerance = 1e-12
    ))
  )
  if (nrow(clone_reference) < 3L * region_size) {
    stop(
      "Only ",
      nrow(clone_reference),
      " clonotypes were found. region_size must be ≤ ",
      floor(nrow(clone_reference) / 3),
      "."
    )
  }
  clone_expression <- obj_sample[[]] |>
    dplyr::filter(
      !is.na(TCR_ID),
      !is.na(CloneID)
    ) |>
    dplyr::group_by(
      TCR_ID,
      CloneID
    ) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(magic_columns),
        ~mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )
  clone_means <- clone_reference |>
    dplyr::left_join(
      clone_expression,
      by = c("TCR_ID", "CloneID")
    )
  rank_tables <- setNames(lapply(magic_columns, function(m) {
    clone_means |>
      dplyr::transmute(
        TCR_ID,
        CloneID,
        CloneFrequency,
        CloneProportion,
        Expression = .data[[m]]
      ) |>
      dplyr::mutate(
        Expression = dplyr::if_else(
          is.finite(Expression) &
            Expression > 0,
          Expression,
          0
        ),
        ZeroExpression = Expression == 0
      ) |>
      dplyr::arrange(
        ZeroExpression,
        dplyr::desc(Expression),
        dplyr::desc(CloneProportion),
        TCR_ID
      ) |>
      dplyr::mutate(
        CloneRankID = dplyr::row_number()
      ) |>
      dplyr::select(-ZeroExpression)
  }), genes)
  region_order <- c(
    "top",
    "middle",
    "low"
  )
  rank_region_tables <- setNames(lapply(
    names(rank_tables),
    function(g) {
      df <- rank_tables[[g]]
      n_clones <- nrow(df)
      middle_start <- floor(
        (n_clones - region_size) / 2
      ) + 1
      middle_end <- middle_start +
        region_size - 1
      dplyr::bind_rows(
        df |>
          dplyr::filter(
            CloneRankID <= region_size
          ) |>
          dplyr::mutate(
            Region = "top"
          ),
        df |>
          dplyr::filter(
            CloneRankID >= middle_start,
            CloneRankID <= middle_end
          ) |>
          dplyr::mutate(
            Region = "middle"
          ),
        df |>
          dplyr::filter(
            CloneRankID > n_clones - region_size
          ) |>
          dplyr::mutate(
            Region = "low"
          )
      )
    }
  ), names(rank_tables))
  rank_region_tcrs <- dplyr::bind_rows(
    rank_region_tables,
    .id = "Gene"
  )
  rank_region_annotations <- rank_region_tcrs |>
    dplyr::group_by(
      Gene,
      Region
    ) |>
    dplyr::summarise(
      xmin = min(CloneRankID),
      xmax = max(CloneRankID),
      x = mean(c(xmin, xmax)),
      .groups = "drop"
    )
  rank_bar_tables <- setNames(lapply(
    names(rank_tables),
    function(g) {
      df <- rank_tables[[g]]
      max_expr <- max(
        df$Expression,
        na.rm = TRUE
      )
      df$ExprScaled <- if (max_expr == 0) {
        0
      } else {
        scales::rescale(
          df$Expression,
          to = c(0, 1),
          from = c(0, max_expr)
        )
      }
      df$GeneColor <- unname(
        gene_colors[g]
      )
      df
    }
  ), names(rank_tables))
  common_ymax <- max(
    clone_reference$CloneProportion,
    na.rm = TRUE
  )
  ybar <- -0.035 * common_ymax
  hbar <- 0.05 * common_ymax
  annotation_offset <- max(
    0.003,
    common_ymax * 0.03
  )
  region_line_y <- max(
    annotation_y,
    common_ymax * 1.05
  )
  region_label_y <- region_line_y +
    annotation_offset
  plot_ymax <- region_label_y +
    annotation_offset * 1.5
  rank_barplots <- setNames(lapply(
    names(rank_bar_tables),
    function(g) {
      df <- rank_bar_tables[[g]]
      region_annotations <- rank_region_annotations |>
        dplyr::filter(Gene == g)
      ggplot(
        df,
        aes(
          CloneRankID,
          CloneProportion
        )
      ) +
        geom_col(
          fill = "black",
          width = 1
        ) +
        geom_tile(
          aes(
            y = ybar,
            fill = ExprScaled
          ),
          height = hbar,
          width = 1
        ) +
        geom_segment(
          data = region_annotations,
          aes(
            x = xmin,
            xend = xmax,
            y = region_line_y,
            yend = region_line_y
          ),
          inherit.aes = FALSE,
          color = "black",
          linewidth = 0.2,
          lineend = "butt"
        ) +
        geom_text(
          data = region_annotations,
          aes(
            x = x,
            y = region_label_y,
            label = Region
          ),
          inherit.aes = FALSE,
          color = "black",
          size = 2.8,
          vjust = 0
        ) +
        scale_fill_gradient(
          low = "white",
          high = unique(df$GeneColor),
          limits = c(0, 1),
          name = "Scaled mean\nexpression"
        ) +
        scale_y_continuous(
          limits = c(
            ybar - hbar / 2,
            plot_ymax
          ),
          expand = c(0, 0)
        ) +
        labs(
          x = "TCR_ID ranked by mean MAGIC-imputed expression",
          y = paste0(
            "CloneProportion within ",
            sample_id_value
          ),
          title = g
        ) +
        coord_cartesian(
          clip = "off"
        ) +
        theme_classic() +
        theme(
          aspect.ratio = 0.25,
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 7),
          axis.ticks.x = element_blank(),
          axis.line = element_line(linewidth = 0.2),
          axis.ticks = element_line(linewidth = 0.2),
          plot.title = element_text(
            size = 7,
            hjust = 0.5
          ),
          legend.key.height = grid::unit(
            0.45,
            "cm"
          ),
          legend.key.width = grid::unit(
            0.35,
            "cm"
          )
        )
    }
  ), names(rank_bar_tables))
  rank_barplot_combined <- patchwork::wrap_plots(
    rank_barplots,
    ncol = barplot_ncol
  )
  marker_order <- names(rank_tables)
  region_top_tcr_by_marker <- rank_region_tcrs |>
    dplyr::group_by(
      Gene,
      Region
    ) |>
    dplyr::arrange(
      CloneRankID,
      .by_group = TRUE
    ) |>
    dplyr::mutate(
      RegionRank = dplyr::row_number()
    ) |>
    dplyr::ungroup()
  region_tcr_overlap <- region_top_tcr_by_marker |>
    dplyr::distinct(
      Gene,
      Region,
      TCR_ID
    ) |>
    dplyr::count(
      Region,
      TCR_ID,
      name = "NumberOfMarkers"
    ) |>
    dplyr::mutate(
      TCRNumber = suppressWarnings(
        as.numeric(sub(
          "^TCR",
          "",
          TCR_ID
        ))
      )
    )
  region_plot_counts <- region_tcr_overlap |>
    dplyr::group_by(Region) |>
    dplyr::summarise(
      PaintedTCRs = sum(
        NumberOfMarkers >=
          minimum_markers_to_color
      ),
      LabelledTCRs = sum(
        NumberOfMarkers >=
          minimum_markers_to_label
      ),
      TotalTCRs = dplyr::n_distinct(TCR_ID),
      .groups = "drop"
    )
  region_top_tcr_by_marker <-
    region_top_tcr_by_marker |>
    dplyr::left_join(
      region_tcr_overlap,
      by = c("Region", "TCR_ID")
    ) |>
    dplyr::mutate(
      PlotRank = dplyr::if_else(
        NumberOfMarkers >=
          minimum_markers_to_color,
        as.numeric(RegionRank),
        NA_real_
      )
    )
  region_tcr_heatmap_tables <- setNames(lapply(
    region_order,
    function(r) {
      region_df <- region_top_tcr_by_marker |>
        dplyr::filter(Region == r)
      region_tcr_order <- region_tcr_overlap |>
        dplyr::filter(Region == r) |>
        dplyr::arrange(
          dplyr::desc(NumberOfMarkers),
          TCRNumber,
          TCR_ID
        )
      tcr_order <- region_tcr_order$TCR_ID
      tidyr::expand_grid(
        TCR_ID = tcr_order,
        Gene = marker_order
      ) |>
        dplyr::left_join(
          region_df |>
            dplyr::select(
              Gene,
              TCR_ID,
              RegionRank,
              NumberOfMarkers,
              PlotRank
            ),
          by = c("TCR_ID", "Gene")
        ) |>
        dplyr::mutate(
          Gene = factor(
            Gene,
            levels = marker_order
          ),
          TCR_ID = factor(
            TCR_ID,
            levels = rev(tcr_order)
          )
        )
    }
  ), region_order)
  gene_positions <- seq_along(
    marker_order
  ) * gene_width
  common_tcr_rows <- max(sapply(
    region_tcr_heatmap_tables,
    function(df) nlevels(df$TCR_ID)
  ))
  region_tcr_heatmaps <- setNames(lapply(
    region_order,
    function(r) {
      df <- region_tcr_heatmap_tables[[r]]
      labelled_tcrs <- region_tcr_overlap |>
        dplyr::filter(
          Region == r,
          NumberOfMarkers >=
            minimum_markers_to_label
        ) |>
        dplyr::pull(TCR_ID)
      region_counts <- region_plot_counts |>
        dplyr::filter(Region == r)
      current_y_levels <- levels(
        df$TCR_ID
      )
      number_of_padding_rows <-
        common_tcr_rows -
        length(current_y_levels)
      padding_y_levels <- if (
        number_of_padding_rows > 0
      ) {
        paste0(
          ".padding_",
          r,
          "_",
          seq_len(number_of_padding_rows)
        )
      } else {
        character(0)
      }
      common_y_levels <- c(
        padding_y_levels,
        current_y_levels
      )
      plot_title <- paste0(
        r,
        "\n(",
        region_counts$LabelledTCRs,
        " top-ranked TCRs in ≥",
        minimum_markers_to_label,
        " markers)"
      )
      ggplot(
        df,
        aes(
          x = as.numeric(Gene) *
            gene_width,
          y = TCR_ID,
          fill = PlotRank
        )
      ) +
        geom_tile(
          width = gene_width,
          height = 0.9,
          color = "grey90",
          linewidth = 0.1
        ) +
        scale_fill_distiller(
          type = "div",
          palette = "RdYlBu",
          direction = 1,
          limits = c(
            1,
            region_size
          ),
          breaks = seq_len(region_size),
          na.value = "white",
          name = "Rank within\nregion",
          guide = guide_colorbar(
            reverse = TRUE
          )
        ) +
        scale_x_continuous(
          breaks = gene_positions,
          labels = marker_order,
          limits = c(
            gene_width / 2,
            (length(marker_order) + 0.5) *
              gene_width
          ),
          expand = c(0, 0)
        ) +
        scale_y_discrete(
          limits = common_y_levels,
          labels = function(x) {
            ifelse(
              x %in% labelled_tcrs,
              x,
              ""
            )
          },
          expand = c(0, 0),
          drop = FALSE
        ) +
        labs(
          x = "Marker",
          y = "TCR_ID",
          title = plot_title
        ) +
        coord_fixed(
          ratio = 1,
          clip = "off"
        ) +
        theme_classic() +
        theme(
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(
            size = 7,
            angle = 90,
            hjust = 1,
            vjust = 0.5
          ),
          axis.text.y = element_text(
            size = 6,
            color = "black"
          ),
          axis.ticks.y = element_blank(),
          axis.line = element_line(linewidth = 0.2),
          axis.ticks = element_line(linewidth = 0.2),
          plot.title = element_text(
            size = 8,
            hjust = 0.5
          ),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7)
        )
    }
  ), region_order)
  region_tcr_heatmap_combined <-
    patchwork::wrap_plots(
      region_tcr_heatmaps,
      ncol = 3,
      widths = rep(
        1,
        length(region_order)
      ),
      guides = "collect"
    ) &
    theme(
      legend.position = "right"
    )
  selected_tcr_region_candidates <- region_tcr_overlap |>
    dplyr::filter(
      NumberOfMarkers >=
        minimum_markers_to_label
    ) |>
    dplyr::mutate(
      RegionPriority = match(
        Region,
        region_order
      )
    )
  region_tcr_conflicts <-
    selected_tcr_region_candidates |>
    dplyr::group_by(TCR_ID) |>
    dplyr::filter(
      dplyr::n_distinct(Region) > 1
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(
      TCR_ID,
      dplyr::desc(NumberOfMarkers),
      RegionPriority
    )
  selected_tcr_regions <-
    selected_tcr_region_candidates |>
    dplyr::arrange(
      TCR_ID,
      dplyr::desc(NumberOfMarkers),
      RegionPriority
    ) |>
    dplyr::distinct(
      TCR_ID,
      .keep_all = TRUE
    ) |>
    dplyr::select(
      TCR_ID,
      Region,
      NumberOfMarkers
    )
  obj_umap <- obj_sub
  obj_umap$TCRRegion <-
    selected_tcr_regions$Region[
      match(
        as.character(obj_umap$TCR_ID),
        selected_tcr_regions$TCR_ID
      )
    ]
  obj_umap$TCRRegion[
    is.na(obj_umap$TCRRegion)
  ] <- "Other"
  obj_umap$TCRRegion <- factor(
    obj_umap$TCRRegion,
    levels = c(
      "Other",
      "low",
      "middle",
      "top"
    )
  )
  selected_umap_coordinates <- as.data.frame(
    SeuratObject::Embeddings(
      obj_umap[["umap"]]
    )
  )
  colnames(selected_umap_coordinates)[1:2] <- c(
    "UMAP_1",
    "UMAP_2"
  )
  selected_umap_coordinates$TCRRegion <-
    obj_umap$TCRRegion[
      match(
        rownames(selected_umap_coordinates),
        colnames(obj_umap)
      )
    ]
  selected_umap_coordinates <-
    selected_umap_coordinates |>
    dplyr::filter(
      TCRRegion != "Other"
    )
  umap_plot <- Seurat::DimPlot(
    obj_umap,
    reduction = "umap",
    group.by = "TCRRegion",
    order = c(
      "Other",
      "low",
      "middle",
      "top"
    ),
    cols = c(
      Other = "grey85",
      top = "#D73027",
      middle = "#FDAE61",
      low = "#4575B4"
    ),
    pt.size = umap_background_size
  ) +
    geom_point(
      data = selected_umap_coordinates,
      aes(
        x = UMAP_1,
        y = UMAP_2,
        color = TCRRegion
      ),
      inherit.aes = FALSE,
      size = umap_selected_size
    ) +
    labs(
      title = paste0(
        sample_id_value,
        ": ",
        sum(selected_tcr_regions$Region == "top"),
        " top TCRs identified in ≥",
        minimum_markers_to_label,
        " markers"
      ),
      color = "TCR region"
    ) +
    theme(
      aspect.ratio = 1,
      axis.line = element_line(
        linewidth = 0.2
      ),
      plot.title = element_text(
        size = 10,
        hjust = 0.5
      )
    )
  compartment_tag <- gsub(
    "[^A-Za-z0-9_-]",
    "_",
    toupper(compartment)
  )
  sample_tag <- gsub(
    "[^A-Za-z0-9_-]",
    "_",
    sample_id_value
  )
  rank_plot_key <- paste0(
    compartment_tag,
    "_",
    sample_tag,
    "_",
    gene_set_name,
    "_rank"
  )
  heatmap_plot_key <- paste0(
    compartment_tag,
    "_",
    sample_tag,
    "_",
    gene_set_name,
    "_heatmap"
  )
  umap_plot_key <- paste0(
    compartment_tag,
    "_",
    sample_tag,
    "_",
    gene_set_name,
    "_umap"
  )
  frequency_key <- paste0(
    compartment_tag,
    "_",
    sample_tag,
    "_",
    gene_set_name
  )
  if (exists(
    "plots_list",
    envir = .GlobalEnv,
    inherits = FALSE
  )) {
    current_plots_list <- get(
      "plots_list",
      envir = .GlobalEnv
    )
    if (!is.list(current_plots_list)) {
      stop("plots_list exists but is not a list.")
    }
  } else {
    current_plots_list <- list()
  }
  current_plots_list[[rank_plot_key]] <-
    rank_barplot_combined
  current_plots_list[[heatmap_plot_key]] <-
    region_tcr_heatmap_combined
  current_plots_list[[umap_plot_key]] <-
    umap_plot
  assign(
    "plots_list",
    current_plots_list,
    envir = .GlobalEnv
  )
  if (exists(
    "freq_lists",
    envir = .GlobalEnv,
    inherits = FALSE
  )) {
    current_freq_lists <- get(
      "freq_lists",
      envir = .GlobalEnv
    )
    if (!is.list(current_freq_lists)) {
      stop("freq_lists exists but is not a list.")
    }
  } else {
    current_freq_lists <- list()
  }
  current_freq_lists[[frequency_key]] <- list(
    sample_id = sample_id_value,
    gene_set = gene_set_name,
    genes = genes,
    region_size = region_size,
    clone_reference = clone_reference,
    clone_means = clone_means,
    rank_tables = rank_tables,
    rank_bar_tables = rank_bar_tables,
    rank_region_tcrs = rank_region_tcrs,
    region_tcr_overlap = region_tcr_overlap,
    region_plot_counts = region_plot_counts,
    region_top_tcr_by_marker =
      region_top_tcr_by_marker,
    region_tcr_heatmap_tables =
      region_tcr_heatmap_tables,
    selected_tcr_regions =
      selected_tcr_regions,
    region_tcr_conflicts =
      region_tcr_conflicts
  )
  assign(
    "freq_lists",
    current_freq_lists,
    envir = .GlobalEnv
  )
  invisible(list(
    sample_id = sample_id_value,
    gene_set = gene_set_name,
    plot_keys = c(
      rank = rank_plot_key,
      heatmap = heatmap_plot_key,
      umap = umap_plot_key
    ),
    frequency_key = frequency_key
  ))
}
phenotype_shifts <- function(data, phenotype, suffix) {
  if (!phenotype %in% colnames(data)) {
    stop("Phenotype column not found: ", phenotype)
  }
  states <- data |>
    dplyr::filter(sample_id %in% annotation_conditions) |>
    dplyr::mutate(state = trimws(as.character(.data[[phenotype]]))) |>
    dplyr::group_by(TCR_ID, sample_id) |>
    dplyr::mutate(total_cells = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(state), state != "") |>
    dplyr::count(TCR_ID, sample_id, total_cells, state, name = "n_state") |>
    dplyr::group_by(TCR_ID, sample_id) |>
    dplyr::filter(n_state == max(n_state)) |>
    dplyr::arrange(state, .by_group = TRUE) |>
    dplyr::summarise(
      annotation = paste(
        sprintf("%s (%d%%)", state, round(100 * n_state / total_cells)),
        collapse = " or "
      ),
      .groups = "drop"
    )
  detected <- data |>
    dplyr::filter(sample_id %in% annotation_conditions) |>
    dplyr::distinct(TCR_ID, sample_id) |>
    dplyr::mutate(detected = TRUE)
  annotations <- tidyr::expand_grid(
    TCR_ID = TCR_priority$TCR_ID,
    sample_id = annotation_conditions
  ) |>
    dplyr::left_join(detected, by = c("TCR_ID", "sample_id")) |>
    dplyr::left_join(states, by = c("TCR_ID", "sample_id")) |>
    dplyr::mutate(
      annotation = dplyr::case_when(
        is.na(detected) ~ "ND",
        is.na(annotation) ~ "Unannotated",
        TRUE ~ annotation
      ),
      column = paste0(sample_id, "_", suffix)
    ) |>
    dplyr::select(TCR_ID, column, annotation) |>
    tidyr::pivot_wider(
      names_from = column,
      values_from = annotation
    )
  annotations
}
negative_phenotype_shifts <- function(data, phenotype, suffix) {
  if (!phenotype %in% colnames(data)) {
    stop("Phenotype column not found: ", phenotype)
  }
  states <- data |>
    dplyr::filter(
      TCR_ID %in% negative_tcr_ids,
      sample_id %in% analysis_conditions
    ) |>
    dplyr::mutate(state = trimws(as.character(.data[[phenotype]]))) |>
    dplyr::group_by(TCR_ID, sample_id) |>
    dplyr::mutate(total_cells = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(state), state != "") |>
    dplyr::count(TCR_ID, sample_id, total_cells, state, name = "n_state") |>
    dplyr::group_by(TCR_ID, sample_id) |>
    dplyr::filter(n_state == max(n_state)) |>
    dplyr::arrange(state, .by_group = TRUE) |>
    dplyr::summarise(
      annotation = paste(
        sprintf("%s (%d%%)", state, round(100 * n_state / total_cells)),
        collapse = " or "
      ),
      .groups = "drop"
    )
  detected <- data |>
    dplyr::filter(
      TCR_ID %in% negative_tcr_ids,
      sample_id %in% analysis_conditions
    ) |>
    dplyr::distinct(TCR_ID, sample_id) |>
    dplyr::mutate(detected = TRUE)
  annotations <- tidyr::expand_grid(
    TCR_ID = negative_tcr_ids,
    sample_id = analysis_conditions
  ) |>
    dplyr::left_join(detected, by = c("TCR_ID", "sample_id")) |>
    dplyr::left_join(states, by = c("TCR_ID", "sample_id")) |>
    dplyr::mutate(
      annotation = dplyr::case_when(
        is.na(detected) ~ "ND",
        is.na(annotation) ~ "Unannotated",
        TRUE ~ annotation
      ),
      column = paste0(sample_id, "_", suffix)
    ) |>
    dplyr::select(TCR_ID, column, annotation) |>
    tidyr::pivot_wider(
      names_from = column,
      values_from = annotation
    )
  annotations
}