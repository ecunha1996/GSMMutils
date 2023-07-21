from GSMMutils.Annotation.genome_annotation import FunctionalAnnotation


def main():
    #genome_annotation = FunctionalAnnotation(r"C:\Users\Bisbii\OneDrive - Universidade do Minho\Algae\Models\Plutheri\blast")
    genome_annotation = FunctionalAnnotation(r"C:\Users\Bisbii\Desktop\synecho")
    genome_annotation.identify_gene_by_homology_from_ec("blastp", "1.14.13.25")



if __name__ == '__main__':
    main()