import unittest

from gsmmutils.annotation import StructuralAnnotation
from os.path import join, dirname

class StructuralAnnotationTest(unittest.TestCase):

    def setUp(self):
        self.gsa = StructuralAnnotation()
        self.data_path = join(dirname(__file__), '../../data/annotation/')

    def test_alignment_evaluation(self):
        genes_annotated, gene_ann_ration = self.gsa.alignment_evaluation(join(self.data_path, 'interproscan_output_test.tsv'), k=3)
        self.assertEquals(gene_ann_ration, 1)
        self.assertEquals(len(genes_annotated), 3)


if __name__ == '__main__':
    unittest.main()
