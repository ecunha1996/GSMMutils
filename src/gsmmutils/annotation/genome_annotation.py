from os import getcwd, chdir, walk
from os.path import getsize
from typing import Union

import pandas as pd

from ..api import UniProt
from ..bio.genome import Genome
from ..utils.utils import run
from Bio.SeqIO import parse


class GenomeAnnotation:
    def __init__(self):
        self._results = {}
        self._genomes = {}
        self._report = {}

    @property
    def genomes(self):
        return self._genomes

    @genomes.setter
    def genomes(self, genomes):
        self._genomes = genomes

    @property
    def results(self):
        return self._results

    @results.setter
    def results(self, results):
        self._results = results

    @property
    def report(self) -> dict:
        return self._report

    @report.setter
    def report(self, report: dict):
        self._report = report

    def load_genes(self, name: str, genome: Union[list, Genome]):
        """
        This function loads the genes.
        Parameters
        ----------
        name: str
            The name of the genome.
        genome: list or Genome
            The list of genes or the genome object.
        Returns
        -------

        """
        self.genomes[name] = genome

    def load_from_fasta(self, filepath: str, name: str = None):
        """
        This function loads the genes from a fasta file.
        Parameters
        ----------
        filepath: str
            The path to the fasta file.

        Returns
        -------

        """
        if not name:
            name = filepath.split('/')[-1].split('.')[0]
        genes = []
        for record in parse(filepath, 'fasta'):
            genes.append(record)
        self.genomes[name] = genes

    def load_genomes_from_folder(self, folder: str):
        """
        This function loads the genomes from a folder.
        Parameters
        ----------
        folder: str
            The path to the folder.

        Returns
        -------

        """
        for dirpath, dirnames, filenames in walk(folder):
            for filename in filenames:
                if filename.endswith('.faa'):
                    genome = Genome()
                    genome.load_from_fasta(f"{dirpath}/{filename}")
                    self.load_genes(filename, genome)

    def load_results(self, path: str, name: str = None, **kwargs):
        """
        This function loads the results.
        Parameters
        ----------
        path: str
            The path to the results.
        name: str
            The name of the results.

        Returns
        -------

        """
        pass

    def calculate_gene_annotation_ratio(self, methods=None):
        """
        This function calculates the gene annotation ratio.
        Returns
        -------

        """
        if methods is None:
            methods = ['interproscan', 'busco']
        for name, result in self.results.items():
            if name not in self.report:
                self.report[name] = {}
            for method in methods:
                genes_with_results = len(result[method].gene_id.unique())
                self.report[name]['Gene Annotation Ratio'] = genes_with_results / len(self.genomes[name])


class StructuralAnnotation(GenomeAnnotation):
    def __init__(self):
        super().__init__()

    def gene_prediction_evaluation(self):
        """
        This function evaluates the gene prediction.
        Returns
        -------

        """
        # TODO: Implement gene prediction evaluation
        pass

    def alignment_evaluation(self, filepath: str = None, k: int = 1):
        """
        This function evaluates the alignment by calculating the gene annotation ratio (genes annotated/total genes).
        Parameters
        ----------
        filepath: str
            The path to the file with the alignment results.
        k: int
            The number of genes in the genome.

        Returns
        -------

        """
        data = pd.read_csv(filepath, sep='\t', header=None)
        unique_genes = len(data[0].unique())
        if not k:
            gene_ann_ration = unique_genes / len(self.genes)
        else:
            gene_ann_ration = unique_genes / k
        return data[0].unique(), gene_ann_ration


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
        uniprot_api = UniProt()
        result = list(uniprot_api.search_by_ec_number(ec_number))
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
