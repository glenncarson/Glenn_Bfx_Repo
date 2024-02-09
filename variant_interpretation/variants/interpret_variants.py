"""This module provides download and variant annotation logic for analyzing clinical variants from VCF files."""

import logging
import os
import subprocess
import pysam


logging.getLogger(__name__)

SNPEFF_PATH = os.getenv("SNPEFF_PATH")
SNPEFF_ANNOTATION_COMMAND_TEMPLATE = "java -Xmx12g -jar {SNPEFF_PATH}SnpSift.jar annotate {annotation_vcf} {vcf_to_annotate} > {output_annotated_vcf}"  # 2> {output_annotated_vcf}.stderr'


class SampleVariants:
    """A class to hold variant information for a sample. It allows variant annotation using a ClinVar annotation file and a method to store vcf entries as records for analysis.
    """
    def __init__(self, sample_id: str, vcf_path: str) -> None:
        self.sample_id = sample_id
        self.sample_vcf = vcf_path
        self.output_annotated_vcf_file_name = None
        self.annotated_records = []

    def annotate_vcf(self, annotations_source_vcf: str, output_path: str, force_overwrite: bool = False) -> None:
        """Annotate the sample with variant annotations from ClinVar
        :param annotations_source_vcf: Path to the ClinVar annotations vcf file.
        :param output_path: The output path for the annotated sample vcf file.
        :param force_overwrite: Whether to overwrite an existing file or not.
        """

        self.output_annotated_vcf_file_name = "".join(
            (output_path, "/", self.sample_id, "_annotated.vcf")
        )
        target_file_exists = os.path.exists(self.output_annotated_vcf_file_name)
        if not target_file_exists or force_overwrite:
            SNPEFF_PATH = os.getenv("SNPEFF_PATH")
            SNPEFF_ANNOTATION_COMMAND_TEMPLATE = "java -Xmx12g -jar {SNPEFF_PATH}SnpSift.jar annotate {annotation_vcf} {vcf_to_annotate} > {output_annotated_vcf}"  # 2> {output_annotated_vcf}.stderr'
            annotation_command = SNPEFF_ANNOTATION_COMMAND_TEMPLATE.format(
                SNPEFF_PATH=SNPEFF_PATH,
                annotation_vcf=annotations_source_vcf,
                vcf_to_annotate=self.sample_vcf,
                output_annotated_vcf=self.output_annotated_vcf_file_name,
            )
            subprocess.run(annotation_command, check=True, shell=True)
        return

    def vcf_to_records(self) -> None:
        """A method to extract records from a vcf file for downstream analysis."""

        with pysam.VariantFile(self.output_annotated_vcf_file_name) as vcf_reader:
            for record in vcf_reader:
                variant_info = {
                    "sample": self.sample_id,
                    "chrom": record.contig,
                    "pos": record.pos,
                    "stop": record.stop,
                    "ref": record.ref,
                    "alts": ", ".join([alt for alt in record.alts]),
                    "qual": record.qual,
                    "filter": ", ".join(filter for filter in record.filter.keys()),
                }
                if "CLNSIG" in record.info.keys():
                    variant_info["clinical_significance"] = record.info["CLNSIG"]
                    variant_info["clinical_variant"] = record.info["CLNVC"]

                self.annotated_records.append(variant_info)
