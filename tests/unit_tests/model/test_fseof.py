import unittest
from src.gsmmutils.model.fseof import FSEOF
from gsmmutils.model import MyModel


class TestFSEOF(unittest.TestCase):

    def setUp(self):
        self.model = MyModel()
        self.fseof = FSEOF(self.model)


if __name__ == '__main__':
    unittest.main()
