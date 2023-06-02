import os
import pathlib
from os import path
from typing import Sequence
from catboost import CatBoostClassifier, cv
import subprocess
import shlex
import pandas as pd

DATA_DIR = path.join(path.curdir, "data")
MODEL_DIR = path.join(path.curdir, "models")


def get_kmers(fasta_path) -> Sequence[str]:
    _ = subprocess.run(
        shlex.split(f'jellyfish count -C -m 31 -s 100M -t 1 {fasta_path}')
    )
    result = subprocess.run(
        shlex.split('jellyfish dump mer_counts.jf'),
        capture_output=True,
    )
    os.remove('mer_counts.jf')

    kmers = set()
    for mer in result.stdout.decode('utf-8').split('\n'):
        if '>' not in mer:
            kmers.add(mer)

    return kmers


def get_matrix_template(antibiotic: str) -> pd.DataFrame:
    return pd.read_csv(path.join(MODEL_DIR, f"{antibiotic}_matrix_head.tsv"), sep="\t")


def fit_matrix(kmers: Sequence[str], template: pd.DataFrame) -> pd.DataFrame:
    x = template
    # Isi satu satu baris baru dengan nol
    x.loc[0] = 0

    for mer in kmers:
        if mer in x.columns:
            x[mer] = 1

    return x


def get_model(antibiotic: str) -> CatBoostClassifier:
    clf = CatBoostClassifier()
    clf.load_model(path.join(MODEL_DIR, f"{antibiotic}_clf.cbm"), format="cbm")
    return clf
