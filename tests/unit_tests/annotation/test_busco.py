import unittest
from GSMMutils.annotation.busco import BuscoDocker
import os
import pytest
IN_GITHUB_ACTIONS = os.getenv("GITHUB_ACTIONS") == "true"

@pytest.mark.skip(reason="Test doesn't work in Github Actions.")
class InterProScanTest(unittest.TestCase):
    def setUp(self):
        self.busco_docker = BuscoDocker(
                                      data_directory="/home/ecunha/GSMMutils/tests/data/annotation",
                                      server="palsson")

    @pytest.mark.skip(reason="Test doesn't work in Github Actions.")
    def test_pull(self):
        self.busco_docker.pull()

    @pytest.mark.skip(reason="Test doesn't work in Github Actions.")
    def test_run(self):
        self.busco_docker.run("-i protein.fasta -o busco_output_test -m proteins -f")


if __name__ == '__main__':
    unittest.main()
