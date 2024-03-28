import unittest
from gsmmutils.annotation.interproscan import InterProScanDocker, InterProScan
import os
import pytest

IN_GITHUB_ACTIONS = os.getenv("GITHUB_ACTIONS") == "true"


@pytest.mark.skip(reason="Test doesn't work in Github Actions.")
class InterProScanTest(unittest.TestCase):
    def setUp(self):
        self.ipsdocker = InterProScanDocker(
            data_directory="/home/ecunha/gsmmutils/tests/data/annotation",
            src_directory="/home/ecunha/gsmmutils/src",
            examples_directory="/home/ecunha/gsmmutils/examples",
            config_directory="/home/ecunha/gsmmutils/config",
            utilities_directory="/home/ecunha/gsmmutils/utilities",
            interproscan_directory="/home/ecunha/interpro/interproscan-5.64-96.0",
            server="palsson")
        self.ips = InterProScan()

    @pytest.mark.skip(reason="Test doesn't work in Github Actions.")
    def test_build(self):
        self.ipsdocker.build()

    @pytest.mark.skip(reason="Test doesn't work in Github Actions.")
    def test_run(self):
        self.ipsdocker.run("-i data/protein.fasta -f TSV -o data/interproscan_output_test.tsv -cpu 8")

    def test_load_results(self):
        self.ips.load_results("data/annotation/interproscan_output_test.tsv")
        self.assertEqual(len(self.ips.results), 1)

    def test_load_results_from_folder(self):
        self.ips.load_results_from_folder("data/annotation")
        self.assertEqual(len(self.ips.results), 1)

    def test_calc_gene_ann_ratio(self):
        self.ips.load_results_from_folder("data/annotation")
        self.ips.load_from_fasta("data/annotation/protein.fasta", "interproscan_output_test")
        self.ips.calculate_gene_annotation_ratio(['interproscan'])
        self.assertEqual(self.ips.report['interproscan_output_test']['Gene Annotation Ratio'], 1.0)


if __name__ == '__main__':
    unittest.main()
