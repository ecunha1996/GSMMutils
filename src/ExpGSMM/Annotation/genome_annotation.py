from os import getcwd, chdir
from os.path import getsize

import pandas as pd
from Bio import SeqIO

from ExpGSMM.api.uniprot import Uniprot
from ExpGSMM.utils.utils import run


class GenomeAnnotation:
    def __init__(self):
        self.genes = None


class StructuralAnnotation(GenomeAnnotation):
    def __init__(self):
        super().__init__()

    def gene_prediction_evaluation(self):
        pass

    def alignment_evaluation(self, method: str, filepath: str = None, k: int=1):
        #data = read_csv(filepath, sep='\t')
        try:
            data = pd.read_csv(filepath, sep='\t', header=None)
        except Exception as e:
            print(e)
            with open(filepath) as f:
                data = f.readlines()
                max_length = max([len(line.split('\t')) for line in data])
                for i in range(len(data)):
                    data[i] = data[i].strip("\n").split('\t')
                    if len(data[i]) < max_length:
                        data[i] += ['-'] * (max_length - len(data[i]))
                data = pd.DataFrame(data)
                data.to_csv(filepath + "_corrected.csv", sep='\t', header=None, index=None)
        unique_genes = len(data[0].unique())
        #gene_counts = data[0].value_counts()
        if not k:
            gene_ann_ration = unique_genes / len(self.genes)
        else:
            gene_ann_ration = unique_genes / k
        print(gene_ann_ration)


class FunctionalAnnotation(GenomeAnnotation):
    def __init__(self, blast_directory):
        super().__init__()
        self.blast_directory = blast_directory


    def identify_gene_by_homology(self, method: str, filepath: str = None):
        """

        :param method:
        :param filepath:
        :return:
        """
        old_dir = getcwd()
        chdir(self.blast_directory)
        if method == 'blastp':
            run(r'blastp -db protein -query "C:\Users\Bisbii\OneDrive - Universidade do Minho\Algae\Models\Dsalina\query.faa" -out '
              r'"C:\Users\Bisbii\OneDrive - Universidade do Minho\Algae\Models\Dsalina\results.txt" -evalue 1  -outfmt 6')
        chdir(old_dir)

    def identify_gene_by_homology_from_ec(self, method: str, ec_number: str):
        """

        :param method:
        :param filepath:
        :return:
        """
        old_dir = getcwd()
        chdir(self.blast_directory)
        # get records from Swissprot, with the ec_number
        uniprot_api = Uniprot()
        result = uniprot_api.search_by_ec_number(ec_number)
        # write the records to a fasta file
        if len(result) > 0:
            with open('query.faa', 'w') as f:
                for record in result:
                    f.write(f">{record['primaryAccession']}\n{record['sequence']['value']}\n")
            if method == 'blastp':
                run(rf'blastp -db protein -query "{self.blast_directory}/query.faa" -out "{self.blast_directory}/results.txt" -evalue 1e-1  -outfmt 6')
            if getsize(f"{self.blast_directory}/results.txt") > 0:
                blast_result = pd.read_csv(f"{self.blast_directory}/results.txt", sep='\t', header=None)
                blast_result.sort_values(by=10, ascending=True, inplace=True)
                print(blast_result.head())
                chdir(old_dir)
            else:
                print(f"No results found with {method}")
        else:
            print("No results found in Swiss-Prot")