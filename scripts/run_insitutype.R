# Check for required packages and install if needed
required_packages <- c("anndataR", "InSituType")
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        if (pkg == "anndataR") {
            pak::pak("scverse/anndataR@v0.1.0")
        } else if (pkg == "InSituType") {
            pak::pak("Nanostring-Biostats/InSituType@v2.0")
        } else {
            install.packages(pkg)
        }
    }
}

suppressPackageStartupMessages({
    library(anndataR)
    library(argparser, quietly = TRUE)
    library(InSituType)
})

parser <- arg_parser("Annotate spatial transcriptomics with insitutype")
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
ref <- read_h5ad(args$ref)

mtx = ref$X
colnames(mtx) = ref$var_names
label = as(ref$obs[[args$ref_column]], "character")
prof = getRNAprofiles(mtx, 0, label)

# use protein if available
if (!is.null(input$obsm) && "intensities" %in% names(input$obsm)) {
    print("using intensities")
    df = input$obsm$intensities

    if ("dummy" %in% colnames(df)) {
        df <- df[, !colnames(df) %in% "dummy"]
    }

    # remove empty columns
    df = df[, colSums(is.na(df)) == 0]

    cohort <- fastCohorting(df, gaussian_transform = TRUE)
} else {
    print("no intensities")
    cohort <- NULL
}

inp_mtx = input$X
colnames(inp_mtx) = input$var_names
inp_mtx = as(inp_mtx, "CsparseMatrix")

res <- insitutypeML(
    x = inp_mtx,
    neg = numeric(nrow(input$X)),
    cohort = cohort,
    reference_profiles = prof
)

names(res$clust) = input$obs_names
write.csv(res$clust, file = args$output)
