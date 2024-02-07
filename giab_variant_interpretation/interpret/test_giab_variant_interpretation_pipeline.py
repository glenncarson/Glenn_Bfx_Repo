import pandas as pd
from pathlib import Path
from giab_variant_interpretation.interpret.interpret_giab_variants import download_variants, download_clinvar_vcf, SampleVariants
# env setup: need to add SNPEFF_PATH as environmental variable, or add snpeff to path
# need to add java runtime env to
args = {
    'output_dir': '/Users/glenn.carson/Coding_Workspace/PycharmProjects/Glenn_Bfx_Repo/analysis',
    'Samples': ['JAS_N36_SNP', 'JAS_P18_SNP', 'M46_SNP', 'M48_SNP'],
    'genome': 'GRCh38',
    'force_overwrite_download_files': False,
    'subdirs' : {'variants_vcf_dir':'raw_vcf_files',
                 'clinvar_annotations_dir': 'clinvar_annotations',
                 'clinvar_annotated_files_dir':'clinvar_annotated_files',
                 'analysis_dir': 'analysis_output'}

}

# make output directories
for dir_name in args['subdirs'].values():
    Path(args['output_dir'], dir_name).mkdir(parents=True, exist_ok=True)

# download variants files. This may take a while.
samples_to_analyze = []
for sample in args['Samples']:
    sample_vcf_path = ''.join((args['output_dir'], '/', args['subdirs']['variants_vcf_dir'], '/', sample, '.vcf')) #removed .gz
    print('Downloading variants for sample {}'.format(sample))
    download_variants(sample, sample_vcf_path, args['force_overwrite_download_files']) #args['genome'],
    samples_to_analyze.append(SampleVariants(sample_id=sample, vcf_path=sample_vcf_path))

# download clinvar file
clinvar_annotations_download_path = ''.join((args['output_dir'], '/', args['subdirs']['clinvar_annotations_dir'], '/{genome}_clinvar.vcf.gz')).format(genome=args['genome'])
download_clinvar_vcf(args['genome'], clinvar_annotations_download_path, force_overwrite=args['force_overwrite_download_files'])

annotated_vcf_path = ''.join((args['output_dir'], '/', args['subdirs']['clinvar_annotated_files_dir']))
all_variant_records = []
for sample in samples_to_analyze:
    sample.annotate_vcf(clinvar_annotations_download_path, annotated_vcf_path, force_overwrite=args['force_overwrite_download_files'])
    sample.vcf_to_records()
    all_variant_records.extend(sample.annotated_records)

all_variants_df = pd.DataFrame.from_records(all_variant_records)
output_df_filename = ''.join((args['output_dir'], '/all_variants.csv'))
all_variants_df.to_csv(output_df_filename, index=False)

# plot 1
