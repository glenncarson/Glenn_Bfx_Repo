import os
import subprocess
import pysam
from giab_variant_interpretation.utils.utils import show_progress
from urllib.request import urlretrieve
from pathlib import Path

ALLOWED_GENOMES = ['GRCh37', 'GRCh38']

# GIAB_URLS = {
#     'NA24385': 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/{genome}/HG002_{genome}_1_22_v4.2.1_benchmark.vcf.gz',
#     'NA24149': 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/{genome}/HG003_{genome}_1_22_v4.2.1_benchmark.vcf.gz',
#     'NA24143': 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/{genome}/HG004_{genome}_1_22_v4.2.1_benchmark.vcf.gz',
#     'NA24631': 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/{genome}/HG005_{genome}_1_22_v4.2.1_benchmark.vcf.gz',
#     'NA24694': 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/ChineseTrio/HG006_NA24694_father/NISTv4.2.1/{genome}/HG006_{genome}_1_22_v4.2.1_benchmark.vcf.gz',
#     'NA24695': 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/ChineseTrio/HG007_NA24695_mother/NISTv4.2.1/{genome}/HG007_{genome}_1_22_v4.2.1_benchmark.vcf.gz',
#     'NA12878': 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/{genome}/HG001_{genome}_1_22_v4.2.1_benchmark.vcf.gz'
# }

# some open source hg38 aligned VCF's: https://figshare.com/articles/dataset/Whole_Exome_Data_VCF_files/13696750
SNP_VCF_URLS = {'JAS_N36_SNP': 'https://figshare.com/ndownloader/files/26347618',
                'JAS_P18_SNP': 'https://figshare.com/ndownloader/files/26347624',
                'M46_SNP': 'https://figshare.com/ndownloader/files/26347630',
                'M48_SNP': 'https://figshare.com/ndownloader/files/26347645'}

SNP_INDEL_URLS = {'JAS_N36_SNP': 'https://figshare.com/ndownloader/files/26347615',
                'JAS_P18_SNP': 'https://figshare.com/ndownloader/files/26347621',
                'M46': 'https://figshare.com/ndownloader/files/26347627',
                'M48': 'https://figshare.com/ndownloader/files/26347633'}

CLINVAR_VCF_URL_TEMPLATE = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_{genome}/clinvar.vcf.gz'

SNPEFF_PATH = os.getenv('SNPEFF_PATH')
SNPEFF_ANNOTATION_COMMAND_TEMPLATE = 'java -Xmx12g -jar {SNPEFF_PATH}SnpSift.jar annotate {annotation_vcf} {vcf_to_annotate} > {output_annotated_vcf}'# 2> {output_annotated_vcf}.stderr'


def download_file_from_url(file_url, output_file_path, force_overwrite=False):
    """Download a file ..."""
    target_file_exists = Path(output_file_path).exists()
    if not target_file_exists or force_overwrite:
        urlretrieve(url=file_url, filename=output_file_path, reporthook=show_progress)
        print('Completed.')
    else:
        print('Skipping download. File already exists.')


def download_vcf_file_and_index_from_url(file_url, output_file_path, get_index_file=False, force_overwrite=False):
    """Download a vcf file and its index from a url."""
    vcf_index_url = ''.join((file_url, ".tbi"))
    vcf_index_output_file_path = ''.join((output_file_path, ".tbi"))
    print(f'Downloading vcf file')
    download_file_from_url(file_url, output_file_path, force_overwrite=force_overwrite)
    if get_index_file:
        print(f'Downloading vcf index file')
        download_file_from_url(vcf_index_url, vcf_index_output_file_path, force_overwrite=force_overwrite)
    return


# add in docstring and types
def download_variants(sample, output_file_path, get_index_file=False, force_overwrite=False): # genome
    # if genome not in ALLOWED_GENOMES:
    #     raise Exception(f'Genome not in allowed list of GIAB genomes for download: {ALLOWED_GENOMES}')
    if sample not in SNP_VCF_URLS:
        raise Exception(f'Sample not in allowed list of samples: {SNP_VCF_URLS.keys()}'
                        f'\nSee: https://figshare.com/articles/dataset/Whole_Exome_Data_VCF_files/13696750')
                        # f'\nSee: https://www.nist.gov/programs-projects/genome-bottle')

    sample_variants_url = SNP_VCF_URLS[sample]#.format(genome=genome)
    download_vcf_file_and_index_from_url(sample_variants_url, output_file_path, get_index_file=get_index_file, force_overwrite=force_overwrite)
    return


def download_clinvar_vcf(genome, output_file_path, force_overwrite=False):
    """Download the clinvar annotations vcf ..."""
    if genome not in ALLOWED_GENOMES:
        raise Exception(f'Genome not in allowed list of GIAB genomes for download: {ALLOWED_GIAB_GENOMES}')
    clinvar_vcf_file_url = CLINVAR_VCF_URL_TEMPLATE.format(genome=genome)
    print(f'Downloading clinvar vcf for {genome}')
    download_vcf_file_and_index_from_url(clinvar_vcf_file_url, output_file_path, force_overwrite)


class SampleVariants:
    def __init__(self, sample_id, vcf_path):
        self.sample_id = sample_id
        self.sample_vcf = vcf_path
        self.output_annotated_vcf_file_name = None
        self.annotated_records = []

    def annotate_vcf(self, annotations_source_vcf, output_path, force_overwrite=False):
        """Annotate the sample with variant annotations from ClinVar"""
        self.output_annotated_vcf_file_name = ''.join((output_path, '/', self.sample_id, '_annotated.vcf'))
        target_file_exists = os.path.exists(self.output_annotated_vcf_file_name)
        if not target_file_exists or force_overwrite:
            SNPEFF_PATH = os.getenv('SNPEFF_PATH')
            SNPEFF_ANNOTATION_COMMAND_TEMPLATE = 'java -Xmx12g -jar {SNPEFF_PATH}SnpSift.jar annotate {annotation_vcf} {vcf_to_annotate} > {output_annotated_vcf}'  # 2> {output_annotated_vcf}.stderr'
            annotation_command = SNPEFF_ANNOTATION_COMMAND_TEMPLATE.format(
                SNPEFF_PATH=SNPEFF_PATH,
                annotation_vcf=annotations_source_vcf,
                vcf_to_annotate=self.sample_vcf,
                output_annotated_vcf=self.output_annotated_vcf_file_name
            )
            subprocess.run(annotation_command, check=True, shell=True)
        return

    def vcf_to_records(self):
        with pysam.VariantFile(self.output_annotated_vcf_file_name) as vcf_reader:
            for record in vcf_reader:
                variant_info = {'sample': self.sample_id,'chrom': record.contig,
                                'pos': record.pos, 'stop': record.stop,
                                'ref': record.ref,'alts': ', '.join([alt for alt in record.alts]),
                                'qual': record.qual, }
                if 'CLNSIG' in record.info.keys():
                    variant_info['clinical_significance'] = record.info['CLNSIG']
                    variant_info['clinical_variant'] = record.info['CLNVC']
                if 'AF_EXAC' in record.info.keys():
                    variant_info['allele_freq_in_EXAC'] = record.info['AF_EXAC']

                self.annotated_records.append(variant_info)
        # vcf_reader = pysam.VariantFile() # or with vcf_reader as


