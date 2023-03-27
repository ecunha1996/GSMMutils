import pandas as pd


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
    def __init__(self):
        super().__init__()