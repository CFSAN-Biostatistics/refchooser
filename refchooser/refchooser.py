# -*- coding: utf-8 -*-

"""This module is part of refchooser.
"""

from __future__ import print_function
from __future__ import absolute_import

from Bio import SeqIO
import gzip
import logging
import os
import pandas as pd
import tempfile

from refchooser import command
from refchooser import utils


def sketch(assemblies, sketch_dir, sketch_size, threads):
    """Create mash sketches to improve the speed of subsequent mash distance calculations.

    Parameters
    ----------
    assemblies : str
        Directory containing assemblies, or a file containing paths to assemblies.
    sketch_dir : str
        Directory where sketches will be stored.
    sketch_size : int
        Each sketch will have at most this many non-redundant min-hashes.
    threads : int
        Number of CPU threads to use.
    """
    if not utils.which("mash"):
        logging.error("Unable to find mash on the path.")
        return

    paths = utils.get_file_list(assemblies)
    if len(paths) == 0:
        return

    utils.mkdir_p(sketch_dir)

    for fasta_path in paths:
        base_file_name = utils.fasta_basename(fasta_path)
        sketch_path = os.path.join(sketch_dir, base_file_name)
        if os.path.isfile(sketch_path + ".msh"):
            logging.info("Skipping already existing %s" % sketch_path + ".msh")
            continue
        command_line = "mash sketch -s %d -p %d -o %s %s" % (sketch_size, threads, sketch_path, fasta_path)
        command.run(command_line)


def get_distance_matrix(sketches):
    """Construct a matrix of mash distances between all pairs of assemblies.

    Parameters
    ----------
    sketches : str
        Directory containing sketches, or a file containing paths to sketches.

    Returns
    -------
    df : Pandas DataFrame
        Matrix of mash distances. Column names and index labels are assembly basenames without extensions.
    """
    if not utils.which("mash"):
        logging.error("Unable to find mash on the path.")
        return

    paths = utils.get_file_list(sketches)
    if len(paths) == 0:
        return

    # Create a file of sketch paths
    with tempfile.NamedTemporaryFile(mode="w") as f_sketches:
        sketch_paths_filename = f_sketches.name
        for sketch_path in paths:
            print(sketch_path, file=f_sketches)
        f_sketches.flush()

        # For each sketch, find the distance to all others
        df_all = pd.DataFrame()
        for sketch_path in paths:
            with tempfile.NamedTemporaryFile() as f_dist:
                dist_filename = f_dist.name
                command_line = "mash dist %s -l %s" % (sketch_path, sketch_paths_filename)
                command.run(command_line, dist_filename)
                df_single = pd.read_csv(dist_filename, sep=None, header=None, usecols=[2], names=["distance"], engine="python")
                assembly_name = utils.basename_no_ext(sketch_path)
                df_all[assembly_name] = df_single["distance"]

        # Change index from integer index to corresponding assembly name
        assembly_index = {i: utils.basename_no_ext(sketch_path) for i, sketch_path in enumerate(paths)}
        df_all.rename(assembly_index, inplace=True)

        return df_all


def distance_matrix(sketches, output_path):
    """Print a matrix of mash distances between all pairs of assemblies and write to a file.

    Parameters
    ----------
    sketches : str
        Directory containing sketches, or a file containing paths to sketches.
    output_path : str
        Path to tab-separated output file.
    """
    if not utils.which("mash"):
        logging.error("Unable to find mash on the path.")
        return

    df = get_distance_matrix(sketches)
    df.to_csv(output_path, sep="\t")


def choose_by_distance(sketches, top_n):
    """Print the list of assemblies having the smallest mash distance to other assemblies
    in a list of assemblies.

    Parameters
    ----------
    sketches : str
        Directory containing sketches, or a file containing paths to sketches.
    top_n : int
        Print the best n candidate references.
    """
    if not utils.which("mash"):
        logging.error("Unable to find mash on the path.")
        return

    df = get_distance_matrix(sketches)
    if len(df) == 0:
        return

    distance_ave = df.mean()
    distance_ave_sorted = distance_ave.sort_values()
    distance_ave_sorted_top = distance_ave_sorted[0: top_n]
    df = pd.DataFrame({"Assembly": distance_ave_sorted_top.index, "Mean_Distance": distance_ave_sorted_top})
    print(df.to_string(index=False))


def choose_by_contigs(assemblies, top_n):
    """Choose an assembly from a collection with minimum number of contigs.

    Parameters
    ----------
    assemblies : str
        Directory containing assemblies, or a file containing paths to assemblies.
    top_n : int
        Print the best n candidate references.
    """
    paths = utils.get_file_list(assemblies)
    if len(paths) == 0:
        return

    # For each assembly, store the number of contigs and file size
    rows = []
    for fasta_path in paths:

        contigs = 0
        size = 0

        if fasta_path.endswith(".gz"):
            open_funct = gzip.open
            mode = "rt"
        else:
            open_funct = open
            mode = "r"
        try:
            with open_funct(fasta_path, mode) as f:
                for seqrecord in SeqIO.parse(f, "fasta"):
                    contigs += 1
                    size += len(seqrecord.seq)
            fasta_base_name = os.path.basename(fasta_path)
            rows.append({"Assembly": fasta_base_name, "Contigs": contigs, "Size": size})
        except FileNotFoundError:
            logging.error("Error opening %s" % fasta_path)

    # Print a list of the best reference genomes sort by number of contigs
    df = pd.DataFrame(rows)
    df.sort_values("Contigs", inplace=True, ascending=True)
    df = df[0: top_n]
    print(df.to_string(index=False))
