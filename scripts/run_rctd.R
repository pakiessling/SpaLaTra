if (!requireNamespace("anndataR", quietly = TRUE)) {
    pak::pak("scverse/anndataR@v0.1.0")
}
if (!requireNamespace("spacexr", quietly = TRUE)) {
    pak::pak("jpromeror/spacexr@33d375cc5d7b7b5db97ee5cc1d3dd32b682afc9e")
}

suppressPackageStartupMessages({
    library(anndataR)
    library(argparser, quietly = TRUE)
    library(spacexr)
    library(SingleCellExperiment)
})

parser <- arg_parser("Annotate spatial transcriptomics with SPLIT")
parser <- add_argument(parser, "--input", help = "Path to h5ad")
parser <- add_argument(parser, "--ref", help = "Path to h5ad")
parser <- add_argument(parser, "--output", help = "Output path for csv")
parser <- add_argument(
    parser,
    "--layer",
    help = "Which matrix to use in the reference",
    default = "X"
)
parser <- add_argument(
    parser,
    "--ref_column",
    help = "Which cell type column to use",
    default = "cell_subtype"
)
parser <- add_argument(
    parser,
    "--max_cores",
    help = "Maximum number of cores for RCTD",
    default = 4L
)
parser <- add_argument(
    parser,
    "--sample_column",
    help = "obs column used to split multi-section h5ads before RCTD (avoids overlapping coordinates)",
    default = NA
)
parser <- add_argument(
    parser,
    "--UMI_min",
    help = "Minimum UMI count per cell passed to create.RCTD (0 = no filtering)",
    default = 0L
)
parser <- add_argument(
    parser,
    "--section",
    help = "Process only cells where sample_column == this value (for parallel mode)",
    default = NA
)

args <- parse_args(parser)

input <- read_h5ad(args$input)
ref <- read_h5ad(args$ref, to = "SingleCellExperiment", x_mapping = "counts")

ref_counts <- counts(ref)
cell_type <- ref[[args$ref_column]]
names(cell_type) <- colnames(ref_counts)

# Replace prohibited characters in cell type names (spacexr disallows '/')
cell_type <- gsub("/", "_", as.character(cell_type), fixed = TRUE)
cell_type <- factor(cell_type)
names(cell_type) <- colnames(ref_counts)

input_counts <- as(input$X, "CsparseMatrix")
input_counts <- Matrix::t(input_counts)
rownames(input_counts) <- input$var_names
colnames(input_counts) <- input$obs_names
coords <- input$obsm$spatial
coords <- as.data.frame(coords)
colnames(coords) <- c("x", "y")
rownames(coords) <- input$obs_names

# remove genes occuring in less than 4 cells to avoid errors in RCTD
ref_counts_filt <- ref_counts[rowSums(ref_counts) > 3, ]

# Determine active cell set (single-section parallel mode vs. full dataset)
if (!is.na(args$section)) {
    if (is.na(args$sample_column)) stop("--section requires --sample_column")
    mask <- as.character(input$obs[[args$sample_column]]) == args$section
    active_counts <- input_counts[, mask, drop = FALSE]
    active_coords <- coords[mask, , drop = FALSE]
} else {
    active_counts <- input_counts
    active_coords <- coords
}

run_rctd_on_cells <- function(
    sub_counts,
    sub_coords,
    ref_counts_filt,
    cell_type,
    max_cores,
    UMI_min = 0
) {
    sub_counts <- sub_counts[rowSums(sub_counts) > 3, , drop = FALSE]
    co_genes <- intersect(rownames(ref_counts_filt), rownames(sub_counts))
    message("[RCTD debug] co_genes: ", length(co_genes))
    message("[RCTD debug] query cells: ", ncol(sub_counts), " | query genes before intersect: ", nrow(sub_counts))
    message("[RCTD debug] ref cells: ", ncol(ref_counts_filt), " | ref genes before intersect: ", nrow(ref_counts_filt))
    sub_ref <- ref_counts_filt[co_genes, ]
    sub_counts <- sub_counts[co_genes, ]
    # Reference() filters cells below min_UMI (default 100) in the co-gene space.
    # With only ~451 co-genes, many cell types lose most of their cells.
    # We set min_UMI = 1 in Reference() and pre-filter to the same threshold so
    # ct_table accurately reflects the cells Reference() will actually keep.
    ref_umi_min <- 1
    ref_cell_totals <- colSums(sub_ref)
    keep_ref_cells <- ref_cell_totals >= ref_umi_min
    n_ref_dropped <- sum(!keep_ref_cells)
    if (n_ref_dropped > 0) {
        message("[RCTD debug] dropping ", n_ref_dropped, " ref cells with 0 UMI in co_genes (",
                ncol(sub_ref) - n_ref_dropped, " remaining)")
        sub_ref   <- sub_ref[, keep_ref_cells, drop = FALSE]
        cell_type <- cell_type[keep_ref_cells]
    }
    ct_table <- table(cell_type)
    rare_types <- names(ct_table[ct_table < 25])
    if (length(rare_types) > 0) {
        message("[RCTD debug] dropping ", length(rare_types), " cell type(s) with < 25 cells: ",
            paste(rare_types, collapse = ", "))
        keep_cells <- !cell_type %in% rare_types
        sub_ref <- sub_ref[, keep_cells, drop = FALSE]
        cell_type <- droplevels(factor(cell_type[keep_cells]))
        ct_table <- table(cell_type)
    }
    message("[RCTD debug] cell types in ref: ", length(ct_table), " | min cells per type: ", min(ct_table), " | max: ", max(ct_table))
    message("[RCTD debug] UMI_min: ", UMI_min, " | counts_MIN: 0 (overridden from default 10)")
    per_cell_counts <- colSums(sub_counts)
    message("[RCTD debug] query counts in co_genes — min: ", min(per_cell_counts),
        " | median: ", median(per_cell_counts), " | max: ", max(per_cell_counts))
    message("[RCTD debug] cells with co_gene counts < 10: ", sum(per_cell_counts < 10))
    reference <- Reference(sub_ref, cell_type, require_int = FALSE, min_UMI = ref_umi_min)
    sample <- SpatialRNA(sub_coords, sub_counts, require_int = FALSE)
    myRCTD <- create.RCTD(
        sample,
        reference,
        max_cores = max_cores,
        UMI_min = UMI_min,
        counts_MIN = 0
    )
    myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
    results <- myRCTD@results$results_df
    message("[RCTD debug] results_df rows: ", if (is.null(results)) "NULL" else nrow(results))
    results
}

# low quality samples tend to crash
tryCatch(
    {
        if (!is.na(args$section)) {
            # single-section parallel mode (dispatched by Snakemake)
            message("[RCTD debug] --- processing section (parallel): ", args$section, " ---")
            message("[RCTD debug] active cells: ", ncol(active_counts))
            result <- run_rctd_on_cells(
                active_counts,
                active_coords,
                ref_counts_filt,
                cell_type,
                args$max_cores,
                args$UMI_min
            )
        } else if (!is.na(args$sample_column)) {
            # original sequential loop (still works if called directly)
            if (!args$sample_column %in% colnames(input$obs)) {
                stop(sprintf(
                    "sample_column '%s' not found in obs. Available columns: %s",
                    args$sample_column,
                    paste(colnames(input$obs), collapse = ", ")
                ))
            }
            sample_ids <- unique(input$obs[[args$sample_column]])
            message("[RCTD debug] sample_column='", args$sample_column, "' | sections found: ", length(sample_ids), " | IDs: ", paste(sample_ids, collapse = ", "))
            parts <- lapply(sample_ids, function(sid) {
                message("[RCTD debug] --- processing section: ", sid, " ---")
                mask <- input$obs[[args$sample_column]] == sid
                run_rctd_on_cells(
                    input_counts[, mask, drop = FALSE],
                    coords[mask, , drop = FALSE],
                    ref_counts_filt,
                    cell_type,
                    args$max_cores,
                    args$UMI_min
                )
            })
            message("[RCTD debug] parts returned: ", length(parts), " | non-null parts: ", sum(!sapply(parts, is.null)))
            result <- do.call(rbind, parts)
        } else {
            result <- run_rctd_on_cells(
                active_counts,
                active_coords,
                ref_counts_filt,
                cell_type,
                args$max_cores,
                args$UMI_min
            )
        }
        if (is.null(result) || nrow(result) == 0) {
            stop("RCTD returned NULL or empty results_df")
        }
        write.csv(result, file = args$output, row.names = TRUE)
    },
    error = function(e) {
        cat("Error occurred:", conditionMessage(e), "\n")
        cat("Saving empty CSV as output due to error.\n")
        error_df <- data.frame(
            spot_class = rep("?", ncol(active_counts)),
            first_type = rep("?", ncol(active_counts)),
            row.names = colnames(active_counts)
        )
        write.csv(error_df, file = args$output, row.names = TRUE)
    }
)
