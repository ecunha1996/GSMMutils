from gsmmutils.annotation.genome_annotation import FunctionalAnnotation


def main():
    genome_annotation = FunctionalAnnotation(r"C:\Users\Bisbii\OneDrive - Universidade do Minho\Algae\Models\Plutheri\blast")
    # genome_annotation = FunctionalAnnotation(r"C:\Users\Bisbii\Desktop\synecho")
    # genome_annotation = FunctionalAnnotation(r"C:\Users\Bisbii\OneDrive - Universidade do Minho\Algae\Models\Ngaditana\blast")
    genome_annotation.identify_gene_by_homology_from_ec("blastp", "2.1.1.17")


# function to order a dictionary by its values
def order_dict_by_values(dictionary):
    """
    This function takes a dictionary as an argument and returns a new dictionary
    that has the same keys and values, but ordered by the values in ascending order.

    Parameters:
    dictionary (dict): The dictionary to be sorted.

    Returns:
    dict: A new dictionary with the same keys and values as the input, but sorted by value.
    """
    return {k: v for k, v in sorted(dictionary.items(), key=lambda item: item[1])}

if __name__ == '__main__':
    main()