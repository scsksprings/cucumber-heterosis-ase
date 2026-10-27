# Data availability

The manuscript identifies NGDC GSA BioProject **PRJCA077972** for 75 paired-end
RNA-seq libraries (150 compressed clean-read files). The project upload log records
transfer completion; archive release/accessibility was not independently verified
during this repository preparation. A GSA CRA accession was not present in the
supplied manuscript or current local submission record and is therefore not invented.

The repository contains STAR gene counts, sample/read filenames, F1 fragment counts,
qPCR Ct observations for the manuscript's 27 genes, phenotype block observations,
annotation mappings, common SNP VCFs and complete statistical outputs.
Raw FASTQs, BAMs and the reference genome are external inputs. The historical
compressed read volume was 227,439,313,849 bytes.

Reference: cucumber Chinese Long v3, with the matching gene annotation used for
STAR, GATK and exon assignment. Gene IDs such as `CsaV3_7G007510` refer to that build.
Use the same annotation snapshot for reproduction; replacing it changes assignment
and enrichment backgrounds.

`docs/source_provenance.json` maps added code and inputs to their original project
relative paths and SHA256 hashes. The original 20261026 archive supplied the figure,
table and manuscript-result files. `checksums.sha256` describes this prepared release.
VCF command-history headers had server-specific filesystem prefixes normalized;
variant/genotype records were unchanged, as recorded in `vcf_header_changes.json`.

The 150 compressed-file MD5 values are in `data/rnaseq/clean_read_checksums.tsv`.
Reference sizes and SHA256 values are in `docs/reference_inputs.tsv`.
