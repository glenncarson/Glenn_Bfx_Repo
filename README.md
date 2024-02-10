# Glenn_Bfx_repo README
## Context
- I made this repository as a very brief public-facing example of the bioinformatics science/engineering I like to do.
- The tools below are meant to be run in two steps:
  - Step One: Download some open-source Exome VCF's and annotate them using ClinVar
  - Step Two: Do some quick analysis. Are the variants good quality? Did we successfully annotate any variants?
- I've attached a pre-run version of my analysis at https://github.com/glenncarson/Glenn_Bfx_Repo/blob/main/analyze_interpreted_variants.md. 
  - Please take a look, but also look at the library code created for the variant-downloading pipeline.


# Tools
## Part One
- `variant_interpretation/bin/download_and_annotate_variants.py`
  - This script download sample vcf files, a clinvar vcf file containing clinical interpretations, annotated the sample vcf records, and then outputs a csv containing all the interpreted variants.
### usage
`python variant_interpretation/bin/download_and_annotate_variants.py /Users/glenn.carson/Coding_Workspace/PycharmProjects/Glenn_Bfx_Repo/analysis --samples JAS_N36 JAS_P18 M46 M48 --genome GRCh38`
`output_path` (required):
- provide an output path to write intermediate and analysis files
`samples` (required): 
- provide these as space-separated names
- options: JAS_N36 JAS_P18 M46 M48
- source of the sample exome SNP vcf's: https://figshare.com/articles/dataset/Whole_Exome_Data_VCF_files/13696750
`genome` (default is GRCh38): 
- Only GRCh38 is allowed at this time, as the sample vcf files are mapped to this build. 
`--force_overwrite_download_files` (default is False):
- if True, overwrites variant files, clinvar annotations, and the annotated vcf files when running the script

## Part Two
- `variant_interpretation/notebooks/analyze_interpreted_variants.ipynb`
- inputs:
  - `variants_csv_path`: path to outputted ready-to-analyze variants from `download_and_annotate_variants.py`
- outputs:
  - Chart one: Assess variant quality scores by sample and group. Check if both groups show statistically different distributions of quality scores.
  - Chart two: Assess the frequency of filtered variants for the 4 samples.
  - Chart three: Assess the variant interpretations for the 4 samples.
  

## Setup
- Clone repository
  - git clone git@github.com:glenncarson/Glenn-Bfx_Repo
- Install python 3.9.18 using pyenv
  - Install pyenv if needed. Instructions: https://github.com/pyenv/pyenv?tab=readme-ov-file
  - set 3.9.18 as your pyenv python
  - pyenv install 3.9.18
  - pyenv global 3.9.18
- find path to python
  - `which python3`
- install poetry
  - `curl -sSL https://install.python-poetry.org | <path-to-python3> - --version 1.7.1`
- navigate to Glenn-Bfx-Repo
  - `cd <Glenn-Bfx-Repo>`
- Create poetry virtual env for the repo
  - `poetry shell`
  - Note: use this command to launch the ve in the future
    - exit with: `exit` command
- Install dependencies
  - `poetry install`
- If using PyCharm, use these instructions to setup. Otherwise, skip to JRE check.
  - Find virtualenv path
    - `which python`
  - navigate to venv dir
    - e.g., `cd /Users/glenn.carson/Coding_Workspace/PycharmProjects/Glenn_Bfx_Repo/.venv`
  - create an excecutable to launch the interpreter
    - `nano python_interp.sh`
    - copy paste these contents:
    - `#!/bin/bash`
    - `source <venv dir>/bin/activate`
      - e.g., `<venv dir>/bin/python "$@"`
    - make the file executable: 
      - `chmod +x python_interp.sh`
    - test the file:
      - `./python_interp.sh`
- Check your Java Runtime Environment (JRE) version
  - This is needed to run SnpEff
  - `java -version`
  - If on an old version, go ahead and update yours. Mine is JRE 17.0.10 at time of writing. 
  - I downloaded my version from here: https://www.oracle.com/java/technologies/javase/jdk17-archive-downloads.html
- Install SnpEff
  - Instructions: https://pcingola.github.io/SnpEff/download/
  - Add it to your poetry venv's path:
  - e.g., /Users/glenn.carson/Coding_Workspace/PycharmProjects/Glenn_Bfx_Repo/.venv/bin/activate
    -  Add this: `export SNPEFF_PATH=/Users/glenn.carson/Coding_Workspace/PycharmProjects/Glenn_Bfx_Repo/snpEff/`
      - 
- Now you are ready to run scripts and notebooks from this repository.