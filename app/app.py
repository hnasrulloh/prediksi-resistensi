import tempfile
import streamlit as st
from os import path
from lib import *

st.title('Prediksi Antibiotik _M. tuberculosis_')

temp_dir = tempfile.mkdtemp()
sample_path = path.join(temp_dir, 'sample.fna')

uploaded = st.file_uploader('Pilih sampel genom')
if uploaded is not None:
    # Simpan sebagai text file
    with open(sample_path, 'w') as f:
        f.write(uploaded.getvalue().decode('utf-8'))

kmers = get_kmers(sample_path)

for antibiotic in ['ethambutol', 'isoniazid', 'rifampin']:
    template = get_matrix_template(antibiotic)
    x = fit_matrix(kmers, template)

    clf = get_model(antibiotic)
    result = 'Resisten' if clf.predict(x)[0] == 1 else 'Sensitif'

    st.markdown(f'Prediksi {antibiotic.title()}: {result}\n' )


# antibiotic = 'isoniazid'
# pheno = 'susceptible'
#
# # kmers = get_kmers(f'examples/H37Rv.fna')
# kmers = get_kmers(f'examples/{antibiotic}_{pheno}.fna')
# template = get_matrix_template(antibiotic)
# x = fit_matrix(kmers, template)
#
# clf = get_model(antibiotic)
# print(clf.predict(x))
