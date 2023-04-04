from ExpAlgae.Annotation.genome_annotation import FunctionalAnnotation


def main():
    genome_annotation = FunctionalAnnotation(r"C:\Users\Bisbii\OneDrive - Universidade do Minho\Algae\Models\Plutheri\blast")
    genome_annotation.identify_gene_by_homology_from_ec("blastp", "1.4.3.16")



if __name__ == '__main__':
    main()