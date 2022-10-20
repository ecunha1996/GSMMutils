import cobra


class Gene(cobra.Gene):
    def __init__(self, gene_id, name=None, chromosome=None, start=None, end=None, strand=None, length=None, description=None):
        self.id = gene_id
        self.name = name
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.strand = strand
        self.length = length
        self.description = description
        super().__init__(gene_id, name)