from gsmmutils.annotation.genome_annotation import StructuralAnnotation


def main():
    ga = StructuralAnnotation()

    ga.alignment_evaluation(r"C:\Users\Bisbii\Downloads\Galaxy3-[InterProScan_on_data_2_(tsv)].tabular", k=10738)
    ga.alignment_evaluation(r"C:\Users\Bisbii\Downloads\Galaxy4-[InterProScan_on_data_1_(tsv)].tabular", k=16587)

    ga.alignment_evaluation(r"C:\Users\Bisbii\Downloads\Galaxy4-[rpsblast_on_data_2].tabular", k=10738)
    ga.alignment_evaluation(r"C:\Users\Bisbii\Downloads\Galaxy5-[rpsblast_on_data_1].tabular", k=16587)


if __name__ == '__main__':
    main()
