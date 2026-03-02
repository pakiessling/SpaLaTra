# Check for required packages and install if needed
if (!requireNamespace("anndataR", quietly = TRUE)) {
    pak::pak("scverse/anndataR@v0.1.0")
}
if (!requireNamespace("PhiSpace", quietly = TRUE)) {
    pak::pak("jiadongm/PhiSpace/pkg@0af720f")
}

suppressPackageStartupMessages({
    library(anndataR)
    library(argparser, quietly = TRUE)
    library(PhiSpace)
    library(scuttle)
})

parser <- arg_parser("Annotate spatial transcriptomics with DOT")
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

args <- parse_args(parser)

input <- read_h5ad(args$input, to = "SingleCellExperiment")
ref <- read_h5ad(args$ref, to = "SingleCellExperiment", x_mapping = "counts")

input <- logTransf(input, assayName = "X", use_log1p = T) # assayName = "counts"
ref <- logTransf(ref, use_log1p = T) # assayName = args$layer,

print("Starting PhiSpace")
PhiRes <- PhiSpaceR_1ref(
    ref, input,
    phenotypes = args$ref_column, PhiSpaceAssay = "data",
    regMethod = "PLS", center = T, scale = F
)
res <- PhiRes$PhiSpaceScore
res <- normPhiScores(res)
write.csv(res, file = args$output, row.names = TRUE)
