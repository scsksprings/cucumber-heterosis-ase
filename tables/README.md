# Publication tables S1–S8

Run from the repository root:

```bash
python tables/build_tables.py
```

The builder reads the processed workbooks in `results/` and the additional ASE
result tables referenced in the code. It writes `tables/output/S1.xlsx` through
`S8.xlsx`, followed by the table index `book.json` and the cell-level presentation
record `formatting_changes.json`.

`headers.py` defines the bilingual column labels. `postprocess_tables.py` applies
the final terminology and rich-text rules after all eight workbooks have been
built. These rules cover explanatory notes and column labels, mathematical
subscripts and italics, the selected gene-symbol labels, and the categorical
qPCR effect labels. For the supplied data, the postprocessing record contains
33 text-cell changes and 649 cells with character-format rules; these sets
overlap.

The postprocessor checks that numeric and Boolean values, formulas, original
identifiers, cell styles, and the remainder of the workbook package are
preserved. S5 has no applicable rule and is preserved byte for byte during this
step. A complete rebuild may change ZIP metadata and unused style indices; the
scientific values and effective styles remain reproducible. Header borders are
set explicitly in the builder to match the supplied tables across pandas
versions.

To apply the same rules to an existing set of eight workbooks:

```bash
python tables/postprocess_tables.py --input-dir path/to/tables --output-dir path/to/final-tables
```

Omit `--output-dir` for an in-place update. Reapplying the rules to the same
formatted workbooks leaves their bytes unchanged. `--manifest` selects a
different path for the cell-level record. Input directories are supplied as
arguments; no machine-specific workspace paths are required.
