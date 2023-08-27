import streamlit as st
from os import path
from lib import *

st.title("Prediksi Antibiotik _M. tuberculosis_ (Global Ed.)")

temp_dir = "temp"
os.makedirs(temp_dir, exist_ok=True)

# This is needed to hold uploaded sample data from user
genome_file = path.join(temp_dir, "genome.fna")

genome_upload = st.file_uploader("Pilih sampel genom")

if genome_upload is not None:
    #
    # Save as text file
    with open(genome_file, "w") as f:
        f.write(genome_upload.getvalue().decode("utf-8"))

    #
    # Decompose
    kmers = decompose_into_kmers(genome_file)

    #
    # Prepare dataset and predict them
    st.markdown("### Hasil Prediksi")

    for antibiotic in ['ethambutol', 'isoniazid', 'rifampin']:
        X = get_table_template(antibiotic)
        X = fit_table(kmers, X)

        clf = get_classifier(antibiotic)
        result = "resisten 🤮" if clf.predict(X)[0] == 1 else "sensitif 🤓"

        st.markdown(f'''
            - **{antibiotic.title()}**  {result.upper()}
        ''')
