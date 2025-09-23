## Running Exomiser and Genomiser

Detailed explanations of the parameters selected and cohort-level run scripts can be found alongside our previous publication:
*An optimized variant prioritization process for rare disease diagnostics: recommendations for Exomiser and Genomiser*
[https://doi.org/10.1186/s13073-025-01546-1]( https://doi.org/10.1186/s13073-025-01546-1)

and its assoicated GitHub repository:
[https://github.com/icooperstein/exomiser_optimization](https://github.com/icooperstein/exomiser_optimization)

### Installation
Installation instructions and instructions to run Exomiser and Genomiser can be found in their documentation: https://exomiser.readthedocs.io/en/latest/running.html
### Filtering VCFs
Before running Exomiser or Genomiser, we recommend applying the following filters to your VCF to remove potential false positive variants:
* GQ ≥ 20
* 0.15 ≤ VAF ≤ 0.85 for heterozygous variants
* ALT != *
* requirements: a "sample.txt" file which simply has the sample_id name (as found in VCF header)

```
bcftools +fill-tags -O z sample.vcf.gz -- -t FORMAT/VAF,HWE | bcftools view -O z -o sample.filtered.vcf.gz -e 'FORMAT/GQ[@sample.txt] < 20 || ( GT[@sample.txt]="het" && ( FORMAT/VAF[@sample.txt] < 0.15 || FORMAT/VAF[@sample.txt] > 0.85  ) ) || ALT ="*"'

```

### Run Exomiser
It is not necessary to replicate our set-up for running Exomiser or Genomiser. \
Necessary files to run Exomiser as a slurm job: 
1. Proband-only or multisample VCF - provided by user
2. Pedigree (for multisample VCF) - provided by user
3. [YAML file](example_yaml.yml) 
4. [application.properties](application.properties)
5. Execution scripts
    - [run_exomiser.sh](run_exomiser.sh)
    - [submit_exomiser.sh](submit_exomiser.sh) \
Linux bash command: ```sbatch submit_exomiser.sh ID run_type``` \
Replace "ID" with your sample ID and "run_type" with naming convention you have named your YML files
