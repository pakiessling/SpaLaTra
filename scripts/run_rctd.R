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
ref_counts <- ref_counts[rowSums(ref_counts) > 3, ]
input_counts <- input_counts[rowSums(input_counts) > 3, ]
co_genes <- intersect(rownames(ref_counts), rownames(input_counts))
ref_counts <- ref_counts[co_genes, ]
input_counts <- input_counts[co_genes, ]


reference <- Reference(ref_counts, cell_type,require_int = FALSE)
sample <- SpatialRNA(coords, input_counts,require_int = FALSE) # support for SPLIT


# low quality samples tend to crash
tryCatch({
  myRCTD <- create.RCTD(sample, reference, max_cores = args$max_cores)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
  rds_path = gsub(".csv", ".rds", args$output)
  saveRDS(myRCTD, file = rds_path)
  result <- myRCTD@results$results_df
  write.csv(result, file = args$output, row.names = TRUE)
}, error = function(e) {

  cat("Error occurred:", conditionMessage(e), "\n")
  cat("Saving empty CSV as output due to error.\n")
  # Save empty CSV
    error_df <- data.frame(
    spot_class = rep("?", ncol(input_counts)),
    first_type = rep("?", ncol(input_counts)),
    row.names = colnames(input_counts)
  )
  write.csv(error_df, file = args$output, row.names = TRUE)
})
