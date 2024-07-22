"""Parse fasta to csv
Clean data to enzymes with available cut site
07/22/2024
Khoa Hoang
"""
import os
from os.path import dirname, basename, isfile, join
from collections import Counter

import matplotlib.pyplot as plt

from SeqUtils import seq_utils

"""
Parse cut site from fasta description header >
"""
def get_cut_site(description: str):
    vals = description.split()
    for v in vals[1:]:
        if set(v).difference(set("ATCG")) == set():
            return v


SCRATCH = os.environ['SCRATCH']
re_data = seq_utils.read_fasta(join(SCRATCH, 'data', 'REBASE', 'All_Type_II_restriction_enzyme_genes_Protein.txt'), "pandas")

re_data["cut_site"] = re_data.Description.apply(get_cut_site)
re_data["cut_site_length"] = re_data["cut_site"].str.len()
re_data["protein_length"] = re_data.Sequence.apply(lambda x: len(x))

re_data = re_data.dropna()
re_data = re_data.drop_duplicates(["Sequence","cut_site"])
cut_site_abundance = re_data["cut_site"].value_counts().to_dict()
re_data["cut_site_abundance"] = re_data.cut_site.map(cut_site_abundance)
re_data["seq"] = [str(i) for i in re_data.Sequence.values]

re_data.to_csv("../data/re.csv", index = None)