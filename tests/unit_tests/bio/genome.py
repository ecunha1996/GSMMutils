import unittest
from gsmmutils.bio.genome import Genome


class TestGenome(unittest.TestCase):
    def setUp(self):
        self.genome = Genome()

    def test_load_from_fasta(self):
        self.genome.from_fasta('data/annotation/protein.fasta')
        self.assertEqual(len(self.genome), 3)


if __name__ == '__main__':
    unittest.main()
