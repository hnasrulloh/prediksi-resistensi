import streamlit as st
from os import path
from lib import *

st.title("Prediksi Antibiotik _M. tuberculosis_")

temp_dir = "temp"
os.makedirs(temp_dir, exist_ok=True)

# This is needed to hold uploaded sample data from user
genome_file = path.join(temp_dir, "genome.fna")
annot_file = path.join(temp_dir, "annot.gff")

genome_upload = st.file_uploader("Pilih sampel genom")
annot_upload = st.file_uploader("Pilih sampel annotasi")

if genome_upload is not None and annot_upload is not None:
    #
    # Save as text file
    with open(genome_file, "w") as f:
        f.write(genome_upload.getvalue().decode("utf-8"))

    with open(annot_file, "w") as f:
        f.write(annot_upload.getvalue().decode("utf-8"))

    #
    # Split genome into genes
    genes_file = os.path.join(temp_dir, "genes.fna")
    split_genome_into_genes(genome_file, annot_file, genes_file)

    #
    # Create Blast DB from genes ref
    ethambutol_genes_ref_file = "gene_refs/ethambutol/ethambutol.fna"
    isoniazid_genes_ref_file = "gene_refs/isoniazid/isoniazid.fna"
    rifampin_genes_ref_file = "gene_refs/rifampin/rifampin.fna"

    db_dir = os.path.join(temp_dir, "db_ref")
    os.makedirs(db_dir, exist_ok=True)

    ethambutol_genes_ref_db = make_blast_db(ethambutol_genes_ref_file, db_dir)
    isoniazid_genes_ref_db = make_blast_db(isoniazid_genes_ref_file, db_dir)
    rifampin_genes_ref_db = make_blast_db(rifampin_genes_ref_file, db_dir)

    #
    # Blastn
    ethambutol_blast_out = f"{genes_file}.ethambutol.blast"
    isoniazid_blast_out = f"{genes_file}.isoniazid.blast"
    rifampin_blast_out = f"{genes_file}.rifampin.blast"

    do_blastn(genes_file, ethambutol_genes_ref_db, ethambutol_blast_out)
    do_blastn(genes_file, isoniazid_genes_ref_db, isoniazid_blast_out)
    do_blastn(genes_file, rifampin_genes_ref_db, rifampin_blast_out)

    #
    # Filter genes
    ethambutol_genes_filtered = get_matches(ethambutol_blast_out, 1e-10)
    isoniazid_genes_filtered = get_matches(isoniazid_blast_out, 1e-10)
    rifampin_genes_filtered = get_matches(rifampin_blast_out, 1e-10)

    ethambutol_filtered_file = f"{genes_file}.ethambutol.filtered"
    isoniazid_filtered_file = f"{genes_file}.isoniazid.filtered"
    rifampin_filtered_file = f"{genes_file}.rifampin.filtered"

    ethambutol_filtered_file = write_seqs(
        ethambutol_genes_filtered, genes_file, desc="ethambutol")
    isoniazid_filtered_file = write_seqs(
        isoniazid_genes_filtered, genes_file, desc="isoniazid")
    rifampin_filtered_file = write_seqs(
        rifampin_genes_filtered, genes_file, desc="rifampin")

    #
    # Decomposition
    ethambutol_kmers = decompose_into_kmers(ethambutol_filtered_file)
    isoniazid_kmers = decompose_into_kmers(isoniazid_filtered_file)
    rifampin_kmers = decompose_into_kmers(rifampin_filtered_file)

    antibiotics_kmers = {
        'ethambutol': ethambutol_kmers,
        'isoniazid': isoniazid_kmers,
        'rifampin': rifampin_kmers,
    }

    #
    # Prepare dataset and predict them
    st.markdown("### Hasil Prediksi")
    for antibiotic in antibiotics_kmers.keys():
        X = get_table_template(antibiotic)
        X = fit_table(antibiotics_kmers[antibiotic], X)

        clf = get_classifier(antibiotic)
        result = "resisten 🤮" if clf.predict(X)[0] == 1 else "sensitif 🤓"

        st.markdown(f'''
            - **{antibiotic.title()}**  {result.upper()}
        ''')
