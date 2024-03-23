import unittest
from unittest.mock import patch, MagicMock
from src.gsmmutils.model.fseof import FSEOF
from gsmmutils import MyModel


class TestFSEOF(unittest.TestCase):

    def setUp(self):
        self.model = MyModel()
        self.fseof = FSEOF(self.model)


if __name__ == '__main__':
    unittest.main()
