from Bio import SeqIO
import pathlib
import uuid


def main(
    input_multifa_dir=snakemake.input[0], OUTPUT=snakemake.params[0],
):
    genome_data = [
        (str(uuid.uuid4()), record)
        for record in SeqIO.parse(input_multifa_dir, format="fasta")
    ]
    # If there is one or less genomes in fasta, exception is raised
    if len(genome_data) <= 1:
        raise Exception(
            f"ERROR: {input_multifa_dir} does not contain more than one sequence"
        )

    # Generates the dir to store the referance
    pathlib.Path(OUTPUT).mkdir(parents=True, exist_ok=True)

    # The first sequence gets a "REF_" prefix to the file name
    for pos, genome in enumerate(genome_data):
        if pos == 0:
            SeqIO.write(genome[1], f"{OUTPUT}/REF_{genome[0]}.fasta", "fasta")
        else:
            SeqIO.write(genome[1], f"{OUTPUT}/{genome[0]}.fasta", "fasta")


if __name__ == "__main__":
    main()
