import unittest
from GSMMutils.annotation.interproscan import InterProScanDocker
IN_GITHUB_ACTIONS = os.getenv("GITHUB_ACTIONS") == "true"



class InterProScanTest(unittest.TestCase):
    def setUp(self):
        self.ipsdocker = InterProScanDocker(
                                      data_directory="/home/ecunha/GSMMutils/tests/data/annotation",
                                      src_directory="/home/ecunha/GSMMutils/src",
                                      examples_directory="/home/ecunha/GSMMutils/examples",
                                      config_directory="/home/ecunha/GSMMutils/config",
                                            utilities_directory="/home/ecunha/GSMMutils/utilities",
                                      interproscan_directory="/home/ecunha/interpro/interproscan-5.64-96.0",
                                    server="palsson")

    @pytest.mark.skipif(IN_GITHUB_ACTIONS, reason="Test doesn't work in Github Actions.")
    def test_build(self):
        self.ipsdocker.build()

    @pytest.mark.skipif(IN_GITHUB_ACTIONS, reason="Test doesn't work in Github Actions.")
    def test_run(self):
        self.ipsdocker.run("-i data/protein.fasta -f TSV -o data/output_test.tsv -cpu 8")


if __name__ == '__main__':
    unittest.main()
