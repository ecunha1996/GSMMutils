from os import getcwd, chdir
from os.path import getsize

import pandas as pd

from GSMMutils.api.uniprot import Uniprot
from GSMMutils.utils.utils import run


class GenomeAnnotation:
    def __init__(self):
        self.genes = None


class StructuralAnnotation(GenomeAnnotation):
    def __init__(self):
        super().__init__()

    def gene_prediction_evaluation(self):
        pass

    def alignment_evaluation(self, filepath: str = None, k: int = 1):
        try:
            data = pd.read_csv(filepath, sep='\t', header=False)
        except Exception as e:
            print(e)
            with open(filepath) as f:
                data_as_list = f.readlines()
                max_length = max([len(line.split('\t')) for line in data_as_list])
                for index, element in enumerate(data_as_list):
                    data_as_list[index] = data_as_list[index].strip("\n").split('\t')
                    if len(data_as_list[index]) < max_length:
                        data_as_list[index] += ['-'] * (max_length - len(data_as_list[index]))
                data_as_list = pd.DataFrame(data_as_list)
                data_as_list.to_csv(f"{filepath}_corrected.csv", sep='\t', header=False, index=False)
                data = data_as_list
        unique_genes = len(data[0].unique())
        if not k:
            gene_ann_ration = unique_genes / len(self.genes)
        else:
            gene_ann_ration = unique_genes / k
        print(gene_ann_ration)


class FunctionalAnnotation(GenomeAnnotation):
    def __init__(self, blast_directory):
        super().__init__()
        self.blast_directory = blast_directory

    def identify_gene_by_homology(self, method: str, query_path: str = None, results_path: str = None):
        """
        This function identifies genes by homology using blastp.
        Parameters
        ----------
        method: str
            The method to be used for the homology search.
        query_path: str
            The path to the query file.
        results_path: str
            The path to the result file.

        Returns
        -------

        """
        old_dir = getcwd()
        chdir(self.blast_directory)
        if method == 'blastp':
            run(rf'blastp -db protein -query {query_path} -out {results_path} -evalue 1 -outfmt 6')
        chdir(old_dir)

    def identify_gene_by_homology_from_ec(self, method: str, ec_number: str):
        """
        This function identifies genes by homology using blastp.
        Parameters
        ----------
        method: str
            The method to be used for the homology search.
        ec_number: str
            The EC number of the enzyme.

        Returns
        -------

        """
        old_dir = getcwd()
        chdir(self.blast_directory)
        # get records from Swiss-prot, with the ec_number
        uniprot_api = Uniprot()
        result = uniprot_api.search_by_ec_number(ec_number)
        # write the records to a fasta file
        if len(result) > 0:
            with open('query.faa', 'w') as f:
                for record in result:
                    f.write(f">{record['primaryAccession']}\n{record['sequence']['value']}\n")
            if method == 'blastp':
                run(rf'blastp -db protein -query "{self.blast_directory}/query.faa" '
                    rf'-out "{self.blast_directory}/results.txt" -evalue 1e-1  -outfmt 6')
            if getsize(f"{self.blast_directory}/results.txt") > 0:
                blast_result = pd.read_csv(f"{self.blast_directory}/results.txt", sep='\t', header=None)
                blast_result.sort_values(by=10, ascending=True, inplace=True)
                print(blast_result.head(20))
                chdir(old_dir)
            else:
                print(f"No results found with {method}")
        else:
            print("No results found in Swiss-Prot")
