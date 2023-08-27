import os
import shlex
import subprocess
import pandas as pd
from typing import List
from catboost import CatBoostClassifier
from Bio import SeqIO


DATA_DIR = os.path.join(os.path.curdir, "data")
MODEL_DIR = os.path.join(os.path.curdir, "models")


def split_genome_into_genes(
    genome_file: str,
    annotation_file: str,
    output_file: str
) -> str:
    """
    Split a genome sequence in FASTA into genes using its RAST annotation
    file from BV-BRC.

    Return output file location.
    """

    # Create a dictionary of genome sequence records
    genome_dict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

    # Open output file for writing
    output_handle = open(output_file, "w")

    # Loop through each line in annotation file
    for line in open(annotation_file):
        # Skip comment lines
        if line.startswith("#"):
            continue

        # Skip empty lines
        if line.strip() == "":
            continue

        # Split line by tabs
        fields = line.strip().split("\t")

        # Extract relevant information
        seqid = fields[0].split('|')[1]  # Adjust to the ID in genome
        feature = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        attributes = fields[8]

        # Check if feature is a gene
        if feature == "CDS":
            # Extract gene ID from attribute string
            [gene_id, *
                gene_desc] = attributes.split("ID=")[1].split("|")[1].split(';')
            # Get the corresponding genome sequence record
            genome_record = genome_dict[seqid]
            # Extract the gene sequence from the genome record
            gene_seq = genome_record.seq[start-1:end]
            # Reverse complement the gene sequence if on minus strand
            if strand == "-":
                gene_seq = gene_seq.reverse_complement()
            # Create a new record for the gene sequence
            gene_record = SeqIO.SeqRecord(
                gene_seq,
                id=gene_id,
                name=gene_id,
                description=";".join(gene_desc)
            )
            # Write the gene record to the output file
            SeqIO.write(gene_record, output_handle, "fasta")

    # Close the output file
    output_handle.close()

    return output_file


def make_blast_db(genes_file, db_dir) -> List[str]:
    """
    Create a Blast DB from a set of genes.

    Return list of DB path. 
    """
    db_path = f"{db_dir}/{os.path.basename(genes_file)}"

    _ = subprocess.run(
        shlex.split(
            f"makeblastdb -in {genes_file} -dbtype nucl -out '{db_path}'")
    )

    return db_path


def do_blastn(genes_file: str, db_path: str, output_file: str) -> str:
    """
    """
    _ = subprocess.run(
        shlex.split(
            f"blastn -task blastn -query {genes_file} -db {db_path} -out {output_file} -outfmt 10")
    )

    return output_file


def get_matches(blast_out: str, evalue_threshold: float) -> pd.DataFrame:
    """
    """
    # Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
    df = pd.read_csv(blast_out, header=None)
    df = df.iloc[:, [0, 1, -2]]
    df.columns = ["query", "ref", "evalue"]
    return df[df.evalue <= evalue_threshold]


def write_seqs(selected_genes: pd.DataFrame, genes_file: str, desc="") -> str:
    """
    """
    seqs = []
    wanted_genes = list(selected_genes['query'])
    for gene in SeqIO.parse(genes_file, "fasta"):
        if gene.name in wanted_genes:
            seqs.append(gene)

    if desc != "":
        desc = f".{desc}"
    selected_genes_file = f"{genes_file}{desc}.filtered"
    SeqIO.write(seqs, selected_genes_file, "fasta")
    return selected_genes_file


def decompose_into_kmers(filtered_genes_file: str) -> List[str]:
    """
    """
    _ = subprocess.run(
        shlex.split(
            f'jellyfish count -C -m 31 -s 100M -t 1 {filtered_genes_file}')
    )
    result = subprocess.run(
        shlex.split('jellyfish dump mer_counts.jf'),
        capture_output=True,
    )
    # os.remove('mer_counts.jf')

    kmers = set()
    for mer in result.stdout.decode('utf-8').split('\n'):
        if '>' not in mer:
            kmers.add(mer)

    return list(kmers)


def get_table_template(antibiotic: str) -> pd.DataFrame:
    """
    """
    return pd.read_csv(os.path.join(MODEL_DIR, f"{antibiotic}_matrix_head.tsv"), sep="\t")


def fit_table(kmers: List[str], template: pd.DataFrame) -> pd.DataFrame:
    """
    Create a single row DataFrame as a new dataset
    """
    X = template
    # Fill the first row with zeros
    X.loc[0] = 0

    # Fill with 1 if the k-mer exists
    for mer in kmers:
        if mer in X.columns:
            X[mer] = 1

    return X


def get_classifier(antibiotic: str) -> CatBoostClassifier:
    """
    """
    clf = CatBoostClassifier()
    clf.load_model(os.path.join(
        MODEL_DIR, f"{antibiotic}_clf.cbm"), format="cbm")
    return clf


# --- TEST DRIVE --- #
#
if __name__ == '__main__':
    # How to split a new genome sample (and its annotation) into a single file of genes.
    os.makedirs("temp/", exist_ok=True)
    genome_file = "examples/83332.12.fna"
    annotation_file = "examples/83332.12.PATRIC.gff"
    genes_file = "temp/genes_83332.12.fna"
    genes_file = split_genome_into_genes(
        genome_file, annotation_file, genes_file)

    # How to make Blast DB from a given set of genes from a new sample.
    genes_ref_file = "gene_refs/ethambutol/ethambutol.fa"
    db_dir = 'temp/db_ref'
    os.makedirs(db_dir, exist_ok=True)
    db_path = make_blast_db(genes_ref_file, db_dir)

    # How to search Blast hit from a given set of genes
    genes_file = genes_file
    db_path = db_path
    blast_out = f"{genes_file}.blast"
    do_blastn(genes_file, db_path, blast_out)

    # How to select genes sequences
    genes_file = genes_file
    filtered_genes = get_matches(blast_out, 1e-10)
    filtered_genes_file = write_seqs(filtered_genes, genes_file)

    # How to decompose genes file into a kmers table
    kmers = decompose_into_kmers(filtered_genes_file)
    X_template = get_table_template("ethambutol")
    X = fit_table(kmers, X_template.copy())

    # How to predict
    clf = get_classifier("ethambutol")
    result = 'Resisten' if clf.predict(X)[0] == 1 else 'Sensitif'
