#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import shutil
from Bio import SeqIO
from .codemltools import model_selection, codeml_settings
import logging


def list_files(work_dir, suffix):
    """Return a list of files with a certain pattern in a directory"""
    file_list = [f for f in os.listdir(work_dir) if f.endswith(suffix)]
    return file_list


def count_dict_values(list_in, dict_in, key):
    """Count nr values from a list having a given key in a dictionary"""
    total_values = 0
    for x in list_in:
        if x in dict_in[key]:
            total_values += 1
    return total_values


def fasta_ids(fasta):
    """Read multi-fasta file and return SPPTAG IDs in a list.
    Format >SPPTAG_geneid"""
    ids = []
    with open(fasta, "r",newline=None) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            ids.append(record.id.split("_")[0])
    return ids


def check_if_duplicates(list_elements):
    """Check if given list contains any duplicates"""
    if len(list_elements) == len(set(list_elements)):
        return False
    else:
        return True


def dup_tags(fasta_file_list):
    """Read a list of fasta files and returns fasta file name
    list containing duplicated genome tags"""
    fasta_dups = []
    for fasta_file in fasta_file_list:
        fasta_handle = SeqIO.parse(fasta_file, "fasta")
        tags = []
        for s in fasta_handle:
            genome_tag = s.id.split("_")[0]
            tags.append(genome_tag)
        # fasta_handle.close()
    if check_if_duplicates(tags):
        fasta_dups.append(fasta_file)
    return fasta_dups


def fasta2phy(msa_input, phy_out):
    """Convert fasta format into format compatible with paml"""
    input_handle = open(msa_input, "r", newline=None)
    output_handle = open(phy_out, "w")
    alignments = SeqIO.parse(input_handle, "fasta")
    msa_seqs = dict()
    for a in alignments:
        seq_id = str(a.id.split("_")[0])
        msa_seqs[seq_id] = str(a.seq)
        bases = len(str(a.seq))
    input_handle.close()
    # write first line of the output file
    output_handle.write("  " + str(len(msa_seqs.keys())) + "  " + str(bases) + "  " + "\n")
    # write sequences in paml format
    for k,v in msa_seqs.items():
        line_out = k + "  " + v + "\n"
        output_handle.write(line_out)
    output_handle.close()


def create_dir(work_dir, dir_name):
    """Create a new folder in the working_dir
    if folder exits, do nothing"""
    dirpath = os.path.join(work_dir, dir_name)
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        pass
        #shutil.rmtree(dirpath)
    else:
        os.mkdir(dirpath)


def move_files(work_dir,file_name,subfolder):
    """Move files to subfolder"""
    file_path = os.path.join(work_dir,file_name)
    shutil.move(file_path,subfolder)


def is_multiple_testing(branches_file):
    """Check if multiple foreground branches"""
    labels = branches_file.keys()
    multiple = 0
    for x in labels:
        if x not in ["0", "1"]:
            multiple += 1
    return multiple


def check_custom(files_path):
    """When custom model selected, checks if .ctl files are in dir
    and returns True if both in dir."""
    st1 = os.path.isfile(os.path.join(files_path, "custom_alt.ctl"))
    st2 = os.path.isfile(os.path.join(files_path, "custom_null.ctl"))
    if st1 and st2:
        return True
    else:
        return False
