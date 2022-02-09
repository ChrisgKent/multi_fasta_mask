from Bio import SeqIO


def main(
    input_multifa_dir=snakemake.input[0],
    OUTPUT=snakemake.output["out_dir"],
    ORDER=snakemake.output["out_order"],
):
    genome_data = [
        (record.id, record) for record in SeqIO.parse(input_multifa_dir, format="fasta")
    ]

    order = []
    for genome in genome_data:
        SeqIO.write(genome[1], f"{OUTPUT}/{genome[0]}.fasta", "fasta")
        order.append(genome[0])

    # Write the order file


if __name__ == "__main__":
    main()
