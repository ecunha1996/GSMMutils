import pandas as pd
from Bio.SeqIO import parse

from ..bio.gene import Gene


class Genome:
    def __init__(self, genes=None, batch=None):
        self._lengths = None
        self._genes = genes or []
        self._genes_ids = [gene.id for gene in self.genes] if genes is not None else []
        if batch is not None:
            self.load(batch)

    @property
    def lengths(self) -> pd.DataFrame:
        return self._lengths

    @lengths.setter
    def lengths(self, value: pd.DataFrame):
        self._lengths = value

    @property
    def genes(self):
        return self._genes

    @genes.setter
    def genes(self, value):
        self._genes = value

    @property
    def genes_ids(self):
        return [gene.id for gene in self.genes]

    @genes_ids.setter
    def genes_ids(self, value):
        self._genes_ids = value

    def __len__(self):
        return len(self.genes)

    def load(self, batch):
        self.genes = [Gene(gene_id=gene_id, length=length) for gene_id, length in batch.to_dict().items()]
        self.lengths = batch

    def set_pathways(self, pathway_map):
        for gene in self.genes:
            if gene.id in pathway_map:
                gene.pathways = pathway_map[gene.id]

    def from_fasta(self, filepath: str):
        """
        This function loads the genes from a fasta file.
        Parameters
        ----------
        filepath: str
            The path to the fasta file.
        """
        with open(filepath) as handle:
            for record in parse(handle, 'fasta'):
                self.genes.append(Gene(gene_id=record.id, sequence=str(record.seq)))

