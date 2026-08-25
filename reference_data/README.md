# Internal SimPheny reference assets

These files are internal SimPheny resources. A normal user does not provide
them at run time.

Expected structure:

```text
reference_data/
├── gene_background.tsv
├── udn/
│   ├── reference.tsv
│   └── phenotype_background.tsv
├── clinvar/reference.tsv
├── phenopacket/reference.tsv
├── decipher/reference.tsv
└── orphanet/reference.tsv
```

Each `reference.tsv` contains:

```text
ID	HPO_Terms	Diagnostic_Genes
```

## Phenotype backgrounds

For **UDN**, the empirical phenotype null uses the broader UDN phenotype
corpus used in the manuscript analysis, rather than only the diagnosed
reference cohort. Store this internal corpus as:

```text
reference_data/udn/phenotype_background.tsv
```

It may use either:

```text
HPO_Terms
HP:0000001;HP:0000002;...
```

(one semicolon-delimited HPO list per participant) or:

```text
HPO_Term
HP:0000001
HP:0000002
...
```

(one annotation per row). Duplicate HPO annotations are intentionally retained.

For **ClinVar, Phenopacket Store, DECIPHER, and Orphanet**, the phenotype
background is generated automatically from all HPO annotations in that
source's `reference.tsv`.

## Gene background

`gene_background.tsv` stores the empirical candidate-gene corpus used for the
gene-p null model. Duplicate gene observations are intentionally retained.

## Restricted resources

UDN and DECIPHER participant-level data must not be committed to the public
repository. Public-source assets should be handled according to the
redistribution terms of each resource.
