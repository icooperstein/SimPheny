# SimPheny

**SimPheny** is a phenotype-first candidate gene prioritization tool for rare disease diagnosis. It compares the phenotypic profile of an undiagnosed patient with diagnosed reference cases and prioritizes candidate genes supported by phenotypically similar patients.

SimPheny accepts one or more query patients, each represented by Human Phenotype Ontology (HPO) terms and a candidate gene list.

## Installation

Python 3.10+ is recommended.

```bash
git clone https://github.com/icooperstein/SimPheny.git
cd SimPheny
pip install -r requirements.txt
```

## Quick start

Create a tab-delimited query file containing:

```text
ID	HPO_Terms	Candidate_Genes
Patient001	HP:0001250;HP:0001263;HP:0000252	RPL13;RYR1;KIF5B
```

Then run:

```bash
python run_simpheny.py \
  --input patient.tsv \
  --reference orphanet
```

By default, SimPheny performs 10,000 empirical iterations per patient-to-patient match and creates a unique timestamped output directory under `results/`.

---

## Query input

The input file must be a tab-delimited TSV with three columns:

| Column | Description |
| --- | --- |
| `ID` | Unique identifier for the query patient |
| `HPO_Terms` | Semicolon-delimited HPO terms |
| `Candidate_Genes` | Semicolon-delimited candidate gene symbols |

For example:

```text
ID	HPO_Terms	Candidate_Genes
Patient001	HP:0001250;HP:0001263;HP:0000252	RPL13;RYR1;KIF5B
```

### Multiple patients

Multiple query patients can be analyzed in a single run by adding one patient per row:

```text
ID	HPO_Terms	Candidate_Genes
Patient001	HP:0001250;HP:0001263;HP:0000252	RPL13;RYR1;KIF5B
Patient002	HP:0004322;HP:0000478;HP:0001249	SET;MTOR;KMT2D
Patient003	HP:0000924;HP:0002758;HP:0002650	SCN2A;CACNA1A;ATP1A3
```

Patient IDs must be unique within the input file. Each patient is analyzed independently and receives an independent candidate-gene ranking.

---

## Reference datasets

Use `--reference` to select the SimPheny reference dataset used for patient matching.

| Option | Reference | Public installation |
| --- | --- | --- |
| `clinvar` | ClinVar-derived reference cases | Available |
| `phenopacket` | Phenopacket Store reference cases | Available |
| `orphanet` | Orphanet disease descriptions | Available |
| `udn` | Undiagnosed Diseases Network | Restricted participant-level data are not distributed with this repository |
| `decipher` | DECIPHER | Restricted participant-level data are not distributed with this repository |

For a standard public installation, use one of the publicly available reference datasets:

```bash
python run_simpheny.py \
  --input patient.tsv \
  --reference orphanet
```

Multiple locally available references can also be selected:

```bash
python run_simpheny.py \
  --input patient.tsv \
  --reference clinvar phenopacket orphanet
```

### UDN and DECIPHER data access

Installation of SimPheny does **not** provide access to restricted UDN or DECIPHER participant-level data. Access to these datasets is governed by their respective data-use and access requirements and is not granted through this repository.

The **SimPheny.iobio** ([simpheny.iobio.io]) web application provides access to SimPheny patient matching using the restricted reference datasets, without distributing the underlying restricted participant-level data.

The command-line implementation supports UDN and DECIPHER for authorized local installations in which the required reference resources are separately available. Support for querying additional hosted reference datasets directly from the command-line implementation is under development.

---

## Command-line options

The general command is:

```bash
python run_simpheny.py \
  --input <query.tsv> \
  --reference <reference>
```

### `--input`

**Required.**

Path to a tab-delimited query file containing `ID`, `HPO_Terms`, and `Candidate_Genes`.

```bash
--input patient.tsv
```

The file may contain either one patient or multiple patients.

### `--reference`

**Required.**

Reference dataset or datasets to use for patient matching.

For a public installation, available reference datasets are:

```text
clinvar
phenopacket
orphanet
```

For example:

```bash
--reference orphanet
```

Multiple locally available references can be supplied:

```bash
--reference clinvar phenopacket orphanet
```

The code also supports `udn` and `decipher` for authorized installations where those restricted reference files are separately available. Installing SimPheny does not provide or grant access to those data.

### `--iterations`

**Optional. Default: `10000`.**

Number of Monte Carlo iterations used to calculate empirical phenotype and gene p-values for each patient-to-patient match.

```bash
--iterations 10000
```

For a faster exploratory run, a smaller number can be used:

```bash
--iterations 1000
```

Because the empirical p-values are estimated by random sampling, fewer iterations provide lower resolution and greater Monte Carlo variability. The manuscript analyses used 10,000 iterations.

### `--seed`

**Optional.**

Sets the random seed used during empirical sampling.

```bash
--seed 42
```

Supplying a seed makes repeated runs with the same inputs and settings reproducible. If no seed is supplied, empirical p-values and SimPheny Scores may vary slightly between runs because of Monte Carlo sampling.

Example:

```bash
python run_simpheny.py \
  --input patient.tsv \
  --reference orphanet \
  --iterations 10000 \
  --seed 42
```

### `--output-dir`

**Optional.**

Specify a custom directory for output files:

```bash
--output-dir my_results
```

If this option is omitted, SimPheny automatically creates a unique timestamped directory under `results/`.

For a single patient:

```text
results/
└── Patient001_orphanet_20260824_173500/
```

For a multi-patient file:

```text
results/
└── patients_3patients_orphanet_20260824_173500/
```

This prevents subsequent analyses from overwriting earlier results.

---

## Example commands

### Single patient

```bash
python run_simpheny.py \
  --input patient.tsv \
  --reference orphanet
```

### Multiple patients

```bash
python run_simpheny.py \
  --input patients.tsv \
  --reference orphanet
```

### Multiple reference datasets

```bash
python run_simpheny.py \
  --input patient.tsv \
  --reference clinvar phenopacket orphanet
```

### Faster exploratory run

```bash
python run_simpheny.py \
  --input patient.tsv \
  --reference orphanet \
  --iterations 1000
```

### Reproducible run

```bash
python run_simpheny.py \
  --input patient.tsv \
  --reference orphanet \
  --iterations 10000 \
  --seed 42
```

### Custom output directory

```bash
python run_simpheny.py \
  --input patient.tsv \
  --reference orphanet \
  --output-dir my_results
```

---

## How SimPheny works

For each query patient and selected reference dataset, SimPheny:

1. compares the query HPO profile with diagnosed reference cases using **PhenoSimJaccard**;
2. identifies reference cases whose diagnostic gene overlaps the query candidate gene list;
3. calculates an empirical phenotype p-value based on the observed phenotypic similarity;
4. calculates an empirical gene p-value based on the frequency of the matched gene in the candidate-gene background;
5. combines phenotype and gene evidence using dataset-specific Empirical Brown's Method parameters;
6. converts the combined p-value to a **SimPheny Score**;
7. aggregates supporting patient matches for each candidate gene; and
8. ranks matched candidate genes for each query patient.

Candidate genes that do not match a diagnostic gene represented in the selected reference dataset do not receive a SimPheny Score.

### Phenotype null distributions

For ClinVar, Phenopacket Store, DECIPHER, and Orphanet, the phenotype null is generated from the pooled HPO annotations in the corresponding diagnosed reference dataset.

For UDN, the phenotype null uses the broader UDN phenotype corpus used in the SimPheny manuscript analysis rather than only the diagnosed UDN reference cohort.

### Candidate-gene background

The empirical gene p-value is calculated using an internal candidate-gene background distribution. Duplicate gene observations are retained because gene frequency in the candidate-gene corpus is part of the empirical null model.

These background resources are internal SimPheny assets and are not supplied by the user at runtime.

---

## SimPheny Scores and confidence tiers

Higher SimPheny Scores indicate stronger combined evidence from phenotypic similarity and gene-level significance.

The manuscript-defined confidence tiers are:

| Confidence | SimPheny Score |
| --- | ---: |
| High | >= 4.5 |
| Medium | >= 2.5 and < 4.5 |
| Low | < 2.5 |

These thresholds were established from benchmarking against diagnosed cases and are intended to help prioritize matches for review.

---

## Gene-level ranking

A candidate gene can be supported by multiple diagnosed reference cases.

For each matched candidate gene, the **SimPheny Gene Score** is calculated from the two highest-scoring supporting patient matches. If only one reference match is available, that match determines the gene score.

Candidate genes are then ranked independently for each query patient.

---

## Output files

Each run produces three files.

### `simpheny_genes.tsv`

Gene-level candidate ranking.

Important fields include:

| Column | Description |
| --- | --- |
| `Query_ID` | Query patient |
| `Candidate_Rank` | Rank within the query patient |
| `Gene_Hit` | Candidate gene supported by reference matches |
| `SimPheny_Gene_Score` | Gene-level SimPheny Score |
| `Num_Reference_Matches` | Number of supporting reference matches |
| `Best_SimPheny_Score` | Highest individual match score |
| `Best_Reference_Source` | Source of the highest-scoring match |
| `Best_Reference_ID` | Reference identifier for the highest-scoring match |
| `Sources` | Reference dataset(s) supporting the gene |

### `simpheny_matches.tsv`

Patient-to-patient match-level results.

Important fields include the query patient, reference source, reference case, matched gene, PhenoSimJaccard score, empirical phenotype p-value, empirical gene p-value, SimPheny Score, and confidence tier.

For restricted datasets, reference identifiers are masked in normal output.

### `run_metadata.json`

Records analysis settings and run information, including:

- query patient IDs;
- selected reference datasets;
- number of empirical iterations;
- matches per query and reference;
- total number of matches;
- runtime; and
- output directory.

---

## Runtime

Runtime depends primarily on:

- number of query patients;
- number and size of selected reference datasets;
- number of candidate-gene matches; and
- number of empirical iterations.

SimPheny reports progress during execution, including the current patient, reference dataset, match number, elapsed time, and estimated remaining time.

For large cohorts, consider using a smaller number of iterations for exploratory analyses before running the final analysis with 10,000 iterations.

---

## Reproducibility

SimPheny uses Monte Carlo sampling to estimate empirical p-values.

For an exactly reproducible run, specify a seed:

```bash
--seed 42
```

Without a seed, repeated analyses can produce slightly different empirical p-values and SimPheny Scores. This is expected and does not affect the deterministic PhenoSimJaccard calculation or identification of candidate-gene overlaps.

---

## Citation

If you use SimPheny in published work, please cite:

> Cooperstein IB, et al. *[SimPheny manuscript citation to be added upon publication]*.

The complete citation will be added following publication.

---

## Contact

For questions, issues, or feature requests, please use the GitHub Issues page.
