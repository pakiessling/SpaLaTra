we assume that both reference and query datasets contain counts

run like this:

nohup snakemake --executor slurm --workflow-profile slurm --conda-frontend conda --config input="input/xenium_hs" output="result/xenium_hs" ref="references/scrna_hs.h5ad" &

input and query need to have counts in .X
layers are not allowed to contain "counts"
obsm has to have coordinate array under "spatial"