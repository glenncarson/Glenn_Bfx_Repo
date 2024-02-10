"""This script is used to download and annotate variant files from John Sayer's Exome Dataset.
See open-source variant files here: https://figshare.com/articles/dataset/Whole_Exome_Data_VCF_files/13696750
"""

import argparse
import logging
import pandas as pd
from pathlib import Path
from variant_interpretation.variants.interpret_variants import SampleVariants
from variant_interpretation.variants.download import (
    download_variants,
    download_clinvar_vcf,
)

logging.basicConfig(level=logging.INFO, format="%(filename)s:%(lineno)d - %(message)s")

SUB_DIRS = {
    "variants_vcf_dir": "raw_vcf_files",
    "clinvar_annotations_dir": "clinvar_annotations",
    "clinvar_annotated_files_dir": "clinvar_annotated_files",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="A script to download variant files and annotate them with ClinVar annotations"
    )
    parser.add_argument("output_path", type=str, help="Path to the output")
    parser.add_argument(
        "-s",
        "--samples",
        nargs="+",
        required=True,
        type=str,
        help="Space-separated samples to annotate. Options: JAS_N36, JAS_P18, M46, M48",
    )
    parser.add_argument(
        "-g",
        "--genome",
        default="GRCh38",
        type=str,
        choices=["GRCh38"],
        help="Genome to annotate. Options: GRCh38",
    )
    parser.add_argument("--force_overwrite_download_files", action="store_true")

    args = parser.parse_args()
    return args


def main(args):
    logging.info("Creating output paths for downloaded variant files")
    for dir_name in SUB_DIRS.values():
        Path(args.output_path, dir_name).mkdir(parents=True, exist_ok=True)

    logging.info("Downloading sample variant files.")
    samples_to_analyze = []
    for sample in args.samples:
        sample_vcf_path = "".join(
            (args.output_path, "/", SUB_DIRS["variants_vcf_dir"], "/", sample, ".vcf")
        )
        download_variants(
            sample, sample_vcf_path, force_overwrite=args.force_overwrite_download_files
        )
        samples_to_analyze.append(
            SampleVariants(sample_id=sample, vcf_path=sample_vcf_path)
        )

    logging.info("Downloading ClinVar annotations.")
    clinvar_annotations_download_path = "".join(
        (
            args.output_path,
            "/",
            SUB_DIRS["clinvar_annotations_dir"],
            "/{genome}_clinvar.vcf.gz",
        )
    ).format(genome=args.genome)
    download_clinvar_vcf(
        args.genome,
        clinvar_annotations_download_path,
        force_overwrite=args.force_overwrite_download_files,
    )

    logging.info("Annotating sample vcf's with Clinvar annotations.")
    annotated_vcf_path = "".join(
        (args.output_path, "/", SUB_DIRS["clinvar_annotated_files_dir"])
    )
    all_variant_records = []
    for sample in samples_to_analyze:
        sample.annotate_vcf(
            clinvar_annotations_download_path,
            annotated_vcf_path,
            force_overwrite=args.force_overwrite_download_files,
        )
        sample.vcf_to_records()
        all_variant_records.extend(sample.annotated_records)

    logging.info("Writing out merged variant records to file.")
    output_df_filename = "".join((args.output_path, "/all_variants.csv"))
    logging.info(f"Output file at: {output_df_filename}")
    all_variants_df = pd.DataFrame.from_records(all_variant_records)
    all_variants_df.to_csv(output_df_filename, index=False)


if __name__ == "__main__":
    main(parse_args())
