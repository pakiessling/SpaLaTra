# Check for required packages and install if needed
required_packages <- c("anndataR")
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        if (pkg == "anndataR") {
            pak::pak("scverse/anndataR@v0.1.0")
        } else {
            install.packages(pkg)
        }
    }
}
suppressPackageStartupMessages({
    library(anndataR)
    library(argparser, quietly = TRUE)
    library(SingleR)
    library(scuttle)
})
parser <- arg_parser("Annotate spatial transcriptomics with DOT")
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
args <- parse_args(parser)

input <- read_h5ad(args$input)
input$layers$counts <- NULL
input <- input$to_SingleCellExperiment(x_mapping = "counts")
ref <- read_h5ad(args$ref, to = "SingleCellExperiment", x_mapping = "counts")

input <- input[, colSums(counts(input)) > 0]
ref <- ref[, colSums(counts(ref)) > 0]

labels <- ref[[args$ref_column]]

input <- logNormCounts(input)
ref <- logNormCounts(ref)

res <- SingleR(
    test = input,
    ref = ref,
    labels = labels,
    num.threads = -1,
    aggr.ref = TRUE
)
write.csv(res, file = args$output, row.names = TRUE)
