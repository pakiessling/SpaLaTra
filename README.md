# SPaLaTra - SPatial Label Transfer

This pipeline performs label transfer from scRNA datasets to spatial datasets with single cell resolution.

A consensus between:
* Insitutype
* RCTD
* TACCO
* SingleR
* PhiSpace

# Starting the pipeline

```bash 
nohup snakemake --executor slurm --workflow-profile slurm --conda-frontend conda --config input="input/xenium_hs" output="result/xenium_hs" ref="references/scrna_hs.h5ad" &
```

# Input

This will run the annotation on all anndata .h5ad in the input directory. Make sure the .obs.index is unique across datasets.
A single cell reference in anndata format is also required.

Query and reference need to have counts in .X

Layers are not allowed to contain the string "counts"

.obsm has to have coordinate x-y numpy array under "spatial"

# Output

The ouput is one .csv per annotation method and one consensus.csv which contains the annotation for every input cell.
