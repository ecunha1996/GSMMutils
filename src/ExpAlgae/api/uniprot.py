import re

import requests
import json
from requests.adapters import HTTPAdapter, Retry
class Uniprot:
    def __init__(self):
        pass



    def search(self):
        pass

    # def search_by_ec_number(self, ec_number: str):
    #     # Send an HTTP GET request to the Swissprot API
    #     url = f'https://rest.uniprot.org/uniprotkb/search?query=(reviewed:true)%20AND%20(ec:{ec_number})&size=500'
    #     response = requests.get(url)
    #
    #     # Parse the JSON response and extract the relevant information
    #     data = json.loads(response.text)
    #     for result in data['results']:
    #         accession = result['primaryAccession']
    #         protein_name = result['proteinDescription']['recommendedName']['fullName']['value']
    #         sequence = result['sequence']['value']
    #         print(accession, protein_name, sequence)
    #     return data

    def get_next_link(self, headers):
        re_next_link = re.compile(r'<(.+)>; rel="next"')
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

    def get_batch(self, batch_url):
        while batch_url:
            response = self.session.get(batch_url)
            response.raise_for_status()
            total = response.headers["x-total-results"]
            yield response, total
            batch_url = self.get_next_link(response.headers)


    def search_by_ec_number(self, ec_number: str):

        retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
        self.session = requests.Session()
        self.session.mount("https://", HTTPAdapter(max_retries=retries))
        url = f'https://rest.uniprot.org/uniprotkb/search?query=(reviewed:true)%20AND%20(ec:{ec_number})&size=500'
        results = []
        for batch, total in self.get_batch(url):
            data = json.loads(batch.text)
            results += data['results']
        return results