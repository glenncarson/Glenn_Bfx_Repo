import logging

logging.getLogger(__name__)
from pathlib import Path
from urllib.request import urlretrieve
from variant_interpretation.utils.utils import show_progress

ALLOWED_GENOMES = ["GRCh37", "GRCh38"]

# some open source Exome data (hg38 aligned VCF's):
# https://figshare.com/articles/dataset/Whole_Exome_Data_VCF_files/13696750
# These are the URL's for SNP VCF's
SNP_VCF_URLS = {
    "JAS_N36": "https://figshare.com/ndownloader/files/26347618",
    "JAS_P18": "https://figshare.com/ndownloader/files/26347624",
    "M46": "https://figshare.com/ndownloader/files/26347630",
    "M48": "https://figshare.com/ndownloader/files/26347645",
}


CLINVAR_VCF_URL_TEMPLATE = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_{genome}/clinvar.vcf.gz"
)


def download_file_from_url(
    file_url: str, output_file_path: str, force_overwrite: bool = False
) -> None:
    """Download a file from a url and save it
    :param file_url: URL to the file to download
    :param output_file_path: Path to the file to save.
    :param force_overwrite: Whether to overwrite an existing file or not.
    """

    target_file_exists = Path(output_file_path).exists()
    if not target_file_exists or force_overwrite:
        urlretrieve(url=file_url, filename=output_file_path, reporthook=show_progress)
        logging.info("Completed.")
    else:
        logging.info("Skipping download. File already exists.")


def download_vcf_file_and_index_from_url(
    file_url: str,
    output_file_path: str,
    get_index_file: bool = False,
    force_overwrite: bool = False,
) -> None:
    """Download a vcf file and its index from a url.
    :param file_url: URL to the file to save.
    :param output_file_path: Path to the file to download.
    :param get_index_file: Whether to download an index file or not for the vcf.
    :param force_overwrite: Whether to overwrite an existing file or not.
    """

    vcf_index_url = "".join((file_url, ".tbi"))
    vcf_index_output_file_path = "".join((output_file_path, ".tbi"))
    download_file_from_url(file_url, output_file_path, force_overwrite=force_overwrite)
    if get_index_file:
        logging.info(f"Downloading vcf index file")
        download_file_from_url(
            vcf_index_url, vcf_index_output_file_path, force_overwrite=force_overwrite
        )
    return


# add in docstring and types
def download_variants(
    sample: str,
    output_file_path: str,
    get_index_file: bool = False,
    force_overwrite: bool = False,
) -> None:
    """Download a variant file and its index for samples in the SNP_VCF_URLS sample set.
    :param sample: Sample name to download.
    :param output_file_path: Path to save the downloaded file.
    :param get_index_file: Whether to download an index file for the vcf.
    :param force_overwrite: Whether to overwrite an existing file or not.
    """

    if sample not in SNP_VCF_URLS:
        raise Exception(
            f"Sample not in allowed list of samples: {SNP_VCF_URLS.keys()}"
            f"\nSee: https://figshare.com/articles/dataset/Whole_Exome_Data_VCF_files/13696750"
        )
    logging.info(f"Downloading vcf for sample {sample}")
    sample_variants_url = SNP_VCF_URLS[sample]
    download_vcf_file_and_index_from_url(
        sample_variants_url,
        output_file_path,
        get_index_file=get_index_file,
        force_overwrite=force_overwrite,
    )
    return


def download_clinvar_vcf(
    genome: str, output_file_path: str, force_overwrite: bool = False
) -> None:
    """Download the ClinVar annotations vcf for a given genome build.
    :param genome: The genome build to use for downloading the ClinVar annotation file.
    :param output_file_path: The target download path of the ClinVar annotation file.
    :param force_overwrite: Whether to overwrite an existing file or not.
    """

    if genome not in ALLOWED_GENOMES:
        raise Exception(
            f"Genome not in allowed list of GIAB genomes for download: {ALLOWED_GENOMES}"
        )
    clinvar_vcf_file_url = CLINVAR_VCF_URL_TEMPLATE.format(genome=genome)
    logging.info(f"Downloading clinvar vcf for {genome}")
    download_vcf_file_and_index_from_url(
        clinvar_vcf_file_url,
        output_file_path,
        get_index_file=True,
        force_overwrite=force_overwrite,
    )
    return
