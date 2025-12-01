We assume that both reference and query datasets contain counts

run like this:

```bash 
nohup snakemake --executor slurm --workflow-profile slurm --conda-frontend conda --config input="input/xenium_hs" output="result/xenium_hs" ref="references/scrna_hs.h5ad" &
``` 

This will run the annotation on all .h5ad in the input directory. Make sure the .obs.index is unique across datasets.

The ouput is one .csv per annotation method and one consensus.csv which contains the annotation for every input cell.

input and query need to have counts in .X

layers are not allowed to contain "counts"

obsm has to have coordinate array under "spatial"
