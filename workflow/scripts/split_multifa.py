from Bio import SeqIO
import pathlib


def main(
    input_multifa_dir=snakemake.input[0], OUTPUT=snakemake.params[0],
):
    genome_data = [
        (record.id, record) for record in SeqIO.parse(input_multifa_dir, format="fasta")
    ]
    # Generates the dir to store the referance
    pathlib.Path(OUTPUT / "referance").mkdir(parents=True, exist_ok=True)

    for pos, genome in enumerate(genome_data):
        if pos == 0:
            SeqIO.write(genome[1], f"{OUTPUT}/REF_{genome[0]}.fasta", "fasta")
        else:
            SeqIO.write(genome[1], f"{OUTPUT}/{genome[0]}.fasta", "fasta")


if __name__ == "__main__":
    main()
