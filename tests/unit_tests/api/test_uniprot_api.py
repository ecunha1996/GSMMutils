import unittest
from GSMMutils.api.uniprot import Uniprot


class MyTestCase(unittest.TestCase):
    def test_uniprot_api_online_reviewed(self):
        ec_number = "2.7.11.1"
        unprot_api = Uniprot()
        result = unprot_api.search_by_ec_number(ec_number)
        self.assertGreater(len(result), 0)

    def test_uniprot_api_online_unreviewed(self):
        ec_number = "2.7.11.1"
        unprot_api = Uniprot()
        result = unprot_api.search_by_ec_number(ec_number, reviewed=False)
        self.assertGreater(len(result), 0)


if __name__ == '__main__':
    unittest.main()
