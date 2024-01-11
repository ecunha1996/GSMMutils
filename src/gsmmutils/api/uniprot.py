import gzip
import json
import re
import shutil
from os import makedirs
from os.path import join, exists
from typing import List
import wget
import requests
from requests.adapters import HTTPAdapter, Retry
from gsmmutils.utils.utils import DATA_PATH


def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


class Uniprot:
    def __init__(self):
        self.data_dir = join(DATA_PATH, join("databases", "uniprot"))
        if not exists(self.data_dir):
            makedirs(self.data_dir)
        self.session = None
        self.db_file = join(self.data_dir, "uniprot_sprot.dat")
        self.db_file_compressed = join(self.data_dir, "uniprot_sprot.dat.gz")

    def search(self):
        """
        Search the UniProt database
        Returns
        -------
        list
        """
        # TODO: Implement search
        pass

    def get_batch(self, batch_url):
        while batch_url:
            response = self.session.get(batch_url)
            response.raise_for_status()
            total = response.headers["x-total-results"]
            yield response, total
            batch_url = get_next_link(response.headers)

    def search_by_ec_number(self, ec_number: str, reviewed: bool = True) -> List[dict]:
        """
        Search for proteins by EC number
        Parameters
        ----------
        ec_number: str
            EC number to search for
        reviewed: bool
            If True, only reviewed proteins are returned

        Returns
        -------
        List[dict]
            A list of proteins matching the EC number
        """
        return self.search_by_ec_number_online(ec_number, reviewed)

    def search_by_ec_number_online(self, ec_number: str, reviewed: bool = True) -> list:
        """
        Search for proteins by EC number using the UniProt REST API
        Parameters
        ----------
        ec_number
        reviewed

        Returns
        -------
        list:
            A list of proteins matching the EC number
        """
        retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
        self.session = requests.Session()
        self.session.mount("https://", HTTPAdapter(max_retries=retries))
        if reviewed:
            url = f'https://rest.uniprot.org/uniprotkb/search?query=(reviewed:true)%20AND%20(ec:{ec_number})&size=500'
        else:
            url = f'https://rest.uniprot.org/uniprotkb/search?query=(ec:{ec_number})&size=500'
        results = []
        for batch, total in self.get_batch(url):
            data = json.loads(batch.text)
            results += data['results']
        return results

    def download_swiss_prot(self):
        """
        Download the SwissProt database
        Returns
        -------

        """
        if not exists(join(self.data_dir, self.db_file_compressed)):
            url = "ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz"
            wget.download(url, join(self.data_dir, self.db_file_compressed))
        if not exists(join(self.data_dir, self.db_file)):
            with gzip.open(join(self.data_dir, self.db_file_compressed), 'rb') as f_in:
                with open(join(self.data_dir, self.db_file), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
