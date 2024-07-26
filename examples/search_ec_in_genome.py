from gsmmutils.annotation.genome_annotation import FunctionalAnnotation


def main():
    genome_annotation = FunctionalAnnotation(r"C:\Users\Bisbii\OneDrive - Universidade do Minho\Algae\Models\Plutheri\blast")
    # genome_annotation = FunctionalAnnotation(r"C:\Users\Bisbii\Desktop\synecho")
    # genome_annotation = FunctionalAnnotation(r"C:\Users\Bisbii\OneDrive - Universidade do Minho\Algae\Models\Ngaditana\blast")
    # genome_annotation = FunctionalAnnotation(r"C:\Users\Bisbii\OneDrive - Universidade do Minho\Algae\Models\Dsalina\blast")
    res = genome_annotation.identify_gene_by_homology_from_ec("blastp", "2.7.1.107")
    print("#"*100)



if __name__ == '__main__':
    main()
