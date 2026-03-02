# Check for required packages and install if needed
required_packages <- c("anndataR", "spacexr")
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        if (pkg == "anndataR") {
            pak::pak("scverse/anndataR@v0.1.0")
        }
        if (pkg == "spacexr") {
            pak::pak("jpromeror/spacexr@33d375cc5d7b7b5db97ee5cc1d3dd32b682afc9e")
        } else {
            install.packages(pkg, repos = "https://cran.uni-muenster.de/")
        }
    }
}

#pak::pak("dmcable/spacexr")

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
parser <- add_argument(parser, "--layer",
    help = "Which matrix to use in the reference",
    default = "X"
)
parser <- add_argument(parser, "--ref_column",
    help = "Which cell type column to use",
    default = "cell_subtype"
)
parser <- add_argument(parser, "--max_cores",
    help = "Maximum number of cores for RCTD",
    default = 4L
)
parser <- add_argument(parser, "--sample_column",
    help = "obs column used to split multi-section h5ads before RCTD (avoids overlapping coordinates)",
    default = NA
)

args <- parse_args(parser)

input <- read_h5ad(args$input)
ref <- read_h5ad(args$ref, to = "SingleCellExperiment", x_mapping = "counts")

ref_counts <- counts(ref)
cell_type <- ref[[args$ref_column]]
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

run_rctd_on_cells <- function(sub_counts, sub_coords, ref_counts_filt, cell_type, max_cores) {
    sub_counts <- sub_counts[rowSums(sub_counts) > 3, , drop = FALSE]
    co_genes <- intersect(rownames(ref_counts_filt), rownames(sub_counts))
    sub_ref <- ref_counts_filt[co_genes, ]
    sub_counts <- sub_counts[co_genes, ]
    reference <- Reference(sub_ref, cell_type, require_int = FALSE)
    sample <- SpatialRNA(sub_coords, sub_counts, require_int = FALSE)
    myRCTD <- create.RCTD(sample, reference, max_cores = max_cores)
    myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
    myRCTD@results$results_df
}

# low quality samples tend to crash
tryCatch({
    if (!is.na(args$sample_column)) {
        sample_ids <- unique(input$obs[[args$sample_column]])
        parts <- lapply(sample_ids, function(sid) {
            mask <- input$obs[[args$sample_column]] == sid
            run_rctd_on_cells(
                input_counts[, mask, drop = FALSE],
                coords[mask, , drop = FALSE],
                ref_counts_filt, cell_type, args$max_cores
            )
        })
        result <- do.call(rbind, parts)
    } else {
        result <- run_rctd_on_cells(input_counts, coords, ref_counts_filt, cell_type, args$max_cores)
    }
    write.csv(result, file = args$output, row.names = TRUE)
}, error = function(e) {
    cat("Error occurred:", conditionMessage(e), "\n")
    cat("Saving empty CSV as output due to error.\n")
    error_df <- data.frame(
        spot_class = rep("?", ncol(input_counts)),
        first_type = rep("?", ncol(input_counts)),
        row.names = colnames(input_counts)
    )
    write.csv(error_df, file = args$output, row.names = TRUE)
})
