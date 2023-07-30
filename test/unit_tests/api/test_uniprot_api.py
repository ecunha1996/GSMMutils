import unittest
from GSMMutils.api.uniprot import Uniprot


class MyTestCase(unittest.TestCase):
    def test_uniprot_api_online(self):
        ec_number = "2.7.11.1"
        unprot_api = Uniprot()
        result = list(unprot_api.search_by_ec_number(ec_number))
        self.assertEqual(len(result), 2633)


if __name__ == '__main__':
    unittest.main()
