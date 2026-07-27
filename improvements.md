# SPaLaTra Improvement Suggestions

## High Priority

### 1. Add input validation in `combine.py`
- `pd.concat([])` on an empty list raises `ValueError` — add existence checks before concatenating
- After joining method outputs, assert no unexpected NaN rows (index mismatches cause silent data corruption)
- Document or implement tie-breaking when `mode()` returns multiple winners

### 2 Output a nice html report
I want to showcase the results of the annotation.
Lets show the annotation and confidence in a UMAP and show some trees / heatmaps that summarize the different annotations.

### 3. Validate globally unique obs indices
The pipeline requires indices to be unique across **all** query files, but nothing enforces this. A collision causes silent row duplication in the final CSV. Add a validation step early in `combine.py`.

### 4. Ideas from feedback

A collaborator gave me some feedback:
I implemented my own plurality consensus voting script that introduced the second label (if available, from methods outputting soft labels) to the voting pool if there's a tie. I didn't touch any of the code for the actual method calls, as I'm assuming defaults are typically used in this context and that normalizations/transformations are done under the hood by the packages themselves (assumed that because of your assert statement that the data are counts/integers). I also didn't implement any QC (dropping low confidence labels or those with a low delta next, etc) as I'm hoping the consensus voting approach is sufficient. With the plurality voting, it seems the annotations are plausible at face value with just 3% drop out (tied vote or voted as unknown), but with majority voting the drop out is higher (~20%). I also stopped using RCTD because it seemed to struggle with a huge fraction of the cells, and it seems to disagree with the other methods a lot (as per Cohen's Kappa and macro F1 score of all pairs of methods). I'll start using those labels to tune some hyperparams like clustering granularity and then manually confirm the biology of the cluster labels with DEGs, but if I'm doing this very wrongly, please do let me know haha.
Please take this into account and make some improvements.
---

## Medium Priority

### 4. Merge `run_rctd.R` and `run_rctd_1c.R`
These two files are nearly identical, differing only in `max_cores`. Parameterize the core count and eliminate the duplicate.

### 5. Surface SLURM CPU count to scripts
RCTD hardcodes `max_cores = 4` but the SLURM config allocates 5 CPUs — one core is always wasted. Pass the allocation via an environment variable or argument, e.g. `$SLURM_CPUS_PER_TASK`.

### 6. Add `log:` and `resources:` directives to the Snakefile
Currently, logs go to stdout/stderr and resources are only in `slurm/config.yaml`. Declaring both in the Snakefile makes it self-documenting and works for non-SLURM runs too:
```python
log: os.path.join(config["output"], "logs", "singler_{sample}.log")
resources:
    mem_mb=120000,
    cpus=20
```

### 7. Fix typo and split the environment file
`enviroment.yml` → `environment.yml`. Consider splitting into a Python env (`env_python.yml`) and an R env (`env_r.yml`) to reduce conflict surface area and installation time.

---

## Lower Priority / Scientific Quality

### 8. Use `workflow.source_path()` for script paths
Replace `"scripts/run_tacco.py"` with `workflow.source_path("scripts/run_tacco.py")` for portability.

### 9. Replace `ref.obs.Donor` with a safe accessor
`run_tacco.py:74` crashes with a `KeyError` if the reference lacks a `Donor` column. Use `ref.obs.get("Donor")` or make it a configurable argument.

### 10. Report annotation quality metrics
The consensus CSV has no confidence column. Consider outputting:
- **Agreement score**: fraction of methods that agreed on the consensus label
- **Ambiguous cells**: flagged when a tie exists
- Per-method columns retained for downstream inspection

### 11. Generate a conda lock file
For exact reproducibility, generate `conda-lock.yml` from `environment.yml` so collaborators get byte-for-byte the same environment.
