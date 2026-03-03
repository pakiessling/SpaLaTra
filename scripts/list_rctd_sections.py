"""
list_rctd_sections.py  –  enumerate tissue sections for parallel RCTD

Reads .obs from an h5ad file (backed mode, no count matrix loaded),
extracts unique values of --sample_column, sanitizes them so they are
safe Snakemake wildcards, and writes a two-column TSV:

    sanitized_name<TAB>original_name

Sanitisation rules
------------------
- Replace every character that is not [A-Za-z0-9_-] with '_'
- Collapse runs of two or more consecutive '_' into a single '_'
- If two original names map to the same sanitised name, append _2, _3, …
  to the duplicates (first occurrence keeps the bare name).
"""

import argparse
import re
import sys

import anndata


def sanitize(name: str) -> str:
    s = re.sub(r"[^A-Za-z0-9_\-]", "_", str(name))
    s = re.sub(r"_+", "_", s)
    return s


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",         required=True, help="Path to query .h5ad")
    parser.add_argument("--sample_column", required=True, help="obs column with section IDs")
    parser.add_argument("--output",        required=True, help="Output TSV path")
    args = parser.parse_args()

    # backed=True skips loading .X into memory
    adata = anndata.read_h5ad(args.input, backed="r")

    if args.sample_column not in adata.obs.columns:
        print(
            f"ERROR: sample_column '{args.sample_column}' not found. "
            f"Available columns: {list(adata.obs.columns)}",
            file=sys.stderr,
        )
        sys.exit(1)

    originals = list(adata.obs[args.sample_column].unique())
    print(f"Found {len(originals)} unique section(s) in '{args.sample_column}'",
          flush=True)

    # Build sanitised names, resolving collisions
    seen: dict[str, int] = {}
    rows: list[tuple[str, str]] = []
    for orig in originals:
        base = sanitize(orig)
        if base not in seen:
            seen[base] = 1
            sane = base
        else:
            seen[base] += 1
            sane = f"{base}_{seen[base]}"
        rows.append((sane, str(orig)))
        print(f"  {orig!r}  →  {sane!r}")

    with open(args.output, "w") as fh:
        for sane, orig in rows:
            fh.write(f"{sane}\t{orig}\n")

    print(f"Wrote {len(rows)} section(s) to {args.output}")


if __name__ == "__main__":
    main()
