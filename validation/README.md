# Author-provided presentation validation

These six JSON records are copied unchanged from the updated author source archive.
They describe the author's Python 3.12.14 / R 4.3.3 rendering run. Fresh integrated
repository checks are in `../docs/validation/`; see `../docs/VALIDATION.md`.

`input_manifest.json` fingerprints the author input files. The prepared repository
normalizes the phenotype CSV to LF and server paths in two VCF headers; those three
packaging differences do not alter statistical values or VCF records. The remaining
97 listed inputs retain the author SHA256 values.
