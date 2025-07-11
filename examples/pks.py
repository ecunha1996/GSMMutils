from Bio import SeqIO


def load_files(filenames):
    results = []
    record_ids = []
    for file in filenames:
        # read with SeqIO
        print(file)
        record = SeqIO.parse(file, "fasta")
        for r in record:
            if r.id not in record_ids:
                record_ids.append(r.id)
                results.append(r)
    return results


def save_files(records, filename):
    with open(filename, "w") as f:
        SeqIO.write(records, f, "fasta")







if __name__ == '__main__':
    records = load_files([r"C:\Users\Bisbii\OneDrive - Universidade do Minho\Algae\PKS\uniprotkb_iterative_polyketide_synthase_2024_04_03_no_duplicates.fasta",
                r"C:\Users\Bisbii\OneDrive - Universidade do Minho\Algae\PKS\uniprotkb_modular_polyketide_synthase_A_2024_04_03_no_duplicates.fasta",
                r"C:\Users\Bisbii\OneDrive - Universidade do Minho\Algae\PKS\uniprotkb_type_II_polyketide_synthase_A_2024_04_03_no_duplicates.fasta",
                r"C:\Users\Bisbii\OneDrive - Universidade do Minho\Algae\PKS\uniprotkb_type_III_polyketide_synthase_2024_04_03_no_duplicates.fasta",
                r"C:\Users\Bisbii\OneDrive - Universidade do Minho\Algae\PKS\query.fasta"])
    save_files(records, r"C:\Users\Bisbii\OneDrive - Universidade do Minho\Algae\PKS\all.fasta")