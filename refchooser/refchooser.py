# -*- coding: utf-8 -*-

"""This module is part of refchooser.
"""

from __future__ import print_function
from __future__ import absolute_import

from Bio import SeqIO
import collections
import gzip
import logging
import numpy as np
import os
import pandas as pd
import sys
import tempfile

from refchooser import command
from refchooser import utils


Metrics = collections.namedtuple("Metrics", ["assembly_name", "n50", "n90", "num_contigs", "total_length"])


def fasta_paths_to_sketch_paths(fasta_paths, sketch_dir, extension_flag=True):
    """Given a list of fasta paths, return the corresponding sketch paths.

    Parameters
    ----------
    fasta_paths : list of str
        List of fasta paths.
    sketch_dir : str
        Directory where sketches will be stored.
    extension : bool, optional, defaults to True
        Controls whether the file extension ".msh" is appended to each of the sketch paths.

    Returns
    -------
    List of sketch paths in the same order at the fasta paths.
    """
    base_file_names = [utils.fasta_basename(fasta_path) for fasta_path in fasta_paths]
    extension = ".msh" if extension_flag else ""
    sketch_paths = [os.path.join(sketch_dir, base_file_name) + extension for base_file_name in base_file_names]
    return sketch_paths


def sketch(assemblies, sketch_dir, sketch_size, threads=1):
    """Create mash sketches to improve the speed of subsequent mash distance calculations.

    Parameters
    ----------
    assemblies : str, or list of str
        Directory containing assemblies, or a file containing paths to assemblies, or a list of fasta paths
    sketch_dir : str
        Directory where sketches will be stored.
    sketch_size : int
        Each sketch will have at most this many non-redundant min-hashes.
    threads : int, optional, defaults to 1
        Number of CPU threads to use.
    """
    if not utils.which("mash"):
        logging.error("Unable to find mash on the path.")
        return

    if isinstance(assemblies, list):
        fasta_paths = assemblies
    else:
        fasta_paths = utils.get_file_list(assemblies)
        if len(fasta_paths) == 0:
            return

    utils.mkdir_p(sketch_dir)

    sketch_paths_with_extension = fasta_paths_to_sketch_paths(fasta_paths, sketch_dir, extension_flag=True)
    sketch_paths_without_extension = fasta_paths_to_sketch_paths(fasta_paths, sketch_dir, extension_flag=False)

    num_already_existing = 0
    for fasta_path, sketch_path_with_extension, sketch_path_without_extension in zip(fasta_paths, sketch_paths_with_extension, sketch_paths_without_extension):
        if os.path.isfile(sketch_path_with_extension):
            num_already_existing += 1
            logging.debug("Skipping already existing %s" % sketch_path_with_extension)
            continue
        command_line = "mash sketch -s %d -p %d -o %s %s" % (sketch_size, threads, sketch_path_without_extension, fasta_path)
        command.run(command_line)

    if num_already_existing > 0:
        logging.debug("Skipped building %d already existing sketches" % num_already_existing)


def get_distance_matrix(fasta_paths, sketch_dir, sketch_size, threads=1):
    """Construct a matrix of mash distances between all pairs of assemblies.

    Parameters
    ----------
    fasta_paths : list of str
        List of paths to fasta files.
    sketch_dir : str
        Directory where sketches will be stored or reused if already existing.
    sketch_size : int
        Each sketch will have at most this many non-redundant min-hashes.
    threads : int, optional, defaults to 1
        Number of CPU threads to use.

    Returns
    -------
    df : Pandas DataFrame
        Matrix of mash distances. Column names and index labels are assembly basenames without extensions.
    """
    if not utils.which("mash"):
        logging.error("Unable to find mash on the path.")
        return

    # Create any missing sketches
    sketch(fasta_paths, sketch_dir, sketch_size, threads=threads)
    sketch_paths_with_extension = fasta_paths_to_sketch_paths(fasta_paths, sketch_dir, extension_flag=True)

    # Create a file of sketch paths.  This will be an input to mash.
    with tempfile.NamedTemporaryFile(mode="w") as f_sketches:
        sketch_paths_filename = f_sketches.name
        for sketch_path in sketch_paths_with_extension:
            print(sketch_path, file=f_sketches)
        f_sketches.flush()

        # Prepare the lower triangle distances
        with tempfile.NamedTemporaryFile() as f_triangle_output:
            triangle_filename = f_triangle_output.name
            command_line = "mash triangle -p %d -l %s" % (threads, sketch_paths_filename)
            command.run(command_line, triangle_filename)

            # Read the triangle, convert to square numpy matrix
            with open(triangle_filename) as f_triangle_input:
                dim = len(fasta_paths)
                a = np.zeros((dim, dim))
                f_triangle_input.readline()  # skip 1st line
                f_triangle_input.readline()  # skip 2nd line
                idx = 1
                for line in f_triangle_input:
                    tokens = line.split()
                    distances = [float(token) for token in tokens[1:]]  # first token is assembly name
                    a[idx, 0: len(distances)] = distances  # partial row
                    a[0: len(distances), idx] = distances  # partial column
                    idx += 1

    # Get a list of short names to use as column headers and row indexes
    assembly_names = [utils.basename_no_ext(sketch_path) for sketch_path in sketch_paths_with_extension]

    # Convert the numpy matrix to Pandas DataFrame
    df = pd.DataFrame(a, index=assembly_names, columns=assembly_names)

    return df


def distance_matrix(assemblies, sketch_dir, sketch_size, output_path, threads=1):
    """Print a matrix of mash distances between all pairs of assemblies and write to a file.

    Parameters
    ----------
    assemblies : str
        Directory containing assemblies, or a file containing paths to assemblies.
    sketch_dir : str
        Directory where sketches will be stored or reused if already existing.
    sketch_size : int
        Each sketch will have at most this many non-redundant min-hashes.
    threads : int, optional, defaults to 1
        Number of CPU threads to use.
    """
    if not utils.which("mash"):
        logging.error("Unable to find mash on the path.")
        return

    fasta_paths = utils.get_file_list(assemblies)
    if len(fasta_paths) == 0:
        return

    df = get_distance_matrix(fasta_paths, sketch_dir, sketch_size, threads)
    df.to_csv(output_path, sep="\t")


def calc_n50(contig_size_list, fraction=0.50):
    """Calculate the N50 (or any other fraction) quality metric of an assembly.

    Parameters
    ----------
    contig_size_list : list of int
        List of sizes of contigs in an assembly.
    fraction : float, optional, defaults to 0.50
        Fraction of the total assembly length covered by contigs at least as
        large as N(fraction).

    Returns
    -------
    n : int
        Contig size such that the collection of all contigs of that length or longer
        contains at least fraction% of the sum of the lengths of all contigs.

    Examples
    --------
    >>> sizes = [539, 552, 596, 714, 763, 820, 908, 929, 1085, 1258, 1356, 1393, 1559, 1562, 2446, 2585, 3500, 4499, 5178, 5241, 6349, 10797, 12402, 18656, 26327, 32855, 32888, 47542, 86790, 93347, 109033, 146629, 146913, 387385, 445110, 535662, 774001, 842337, 1161370]
    >>> calc_n50(sizes)
    774001
    >>> calc_n50(sizes, 0.50)
    774001
    >>> calc_n50(sizes, 0.90)
    109033
    """
    contig_size_list = np.array(sorted(contig_size_list, reverse=True))
    total_length = np.sum(contig_size_list)
    cumsum = np.cumsum(contig_size_list)
    covered_length = int(total_length * fraction)
    cumsum_long_enough_flags = cumsum >= covered_length
    n = contig_size_list[cumsum_long_enough_flags][0]
    return n


def compute_metrics(fasta_path):
    """Compute metrics for a specified fasta file.

    Parameters
    ----------
    fasta_path : str
        Path to fasta file.

    Returns
    -------
    metrics : Metrics
        Metrics for this file
    """
    contig_size_list = []

    if fasta_path.endswith(".gz"):
        open_funct = gzip.open
        mode = "rt"
    else:
        open_funct = open
        mode = "r"
    try:
        with open_funct(fasta_path, mode) as f:
            for seqrecord in SeqIO.parse(f, "fasta"):
                contig_size_list.append(len(seqrecord.seq))
    except FileNotFoundError:
        logging.error("Error opening %s" % fasta_path)
        return None

    assembly_name = utils.fasta_basename(fasta_path)
    num_contigs = len(contig_size_list)
    total_length = sum(contig_size_list)

    # Compute N50
    n50 = calc_n50(contig_size_list, 0.50)
    n90 = calc_n50(contig_size_list, 0.90)
    return Metrics(assembly_name, n50, n90, num_contigs, total_length)


def metrics(assemblies, sketch_dir, sketch_size, top_n, sort_by=None, threads=1):
    """Print a table of metrics to help choose an assembly from a collection.

    Parameters
    ----------
    assemblies : str
        Directory containing assemblies, or a file containing paths to assemblies.
    sketch_dir : str
        Directory where sketches will be stored or reused if already existing.
    sketch_size : int
        Each sketch will have at most this many non-redundant min-hashes.
    sort_by : str, optional, defaults to None
        Sort results by table column.
    top_n : int
        Print the best n candidate references.
    threads : int, optional, defaults to 1
        Number of CPU threads to use.
    """
    available_sort_columns = ["Assembly", "N50", "N90", "Contigs", "Length", "Mean_Distance", "Path", "Score"]
    if sort_by is None:
        sort_by = "Score"
    sort_by_lower = sort_by.lower()
    sort_by = None
    for column in available_sort_columns:
        if sort_by_lower == column.lower():
            sort_by = column
            break
    if sort_by is None:
        logging.warning("Unrecognized sort column.  Choose from: %s." % ", ".join(available_sort_columns))
        return

    if not utils.which("mash"):
        logging.error("Unable to find mash on the path.")
        return

    fasta_paths = utils.get_file_list(assemblies)
    if len(fasta_paths) == 0:
        return

    # For each assembly, gather its metrics
    fasta_paths_having_metrics = []
    metrics_rows = []
    for fasta_path in fasta_paths:
        metrics = compute_metrics(fasta_path)
        if metrics:
            fasta_paths_having_metrics.append(fasta_path)
            metrics_rows.append(metrics)

    if len(fasta_paths_having_metrics) == 0:
        return

    # For each assembly, compute the average distance to all others
    df = get_distance_matrix(fasta_paths_having_metrics, sketch_dir, sketch_size, threads)
    mean_distance_series = df.sum() / (len(df) - 1)  # Don't include zero distance to self in denominator

    # Build the composite dataframe joining the metrics and distance calculations
    n50_list = [metrics.n50 for metrics in metrics_rows]
    n90_list = [metrics.n90 for metrics in metrics_rows]
    num_contigs_list = [metrics.num_contigs for metrics in metrics_rows]
    total_length_list = [metrics.total_length for metrics in metrics_rows]
    df = pd.DataFrame({"Assembly": mean_distance_series.index, "N50": n50_list, "N90": n90_list, "Contigs": num_contigs_list, "Length": total_length_list, "Mean_Distance": mean_distance_series, "Path": fasta_paths_having_metrics})

    # Compute a "score"
    df["Score"] = df.N50 / df.Mean_Distance

    # Force the order of the columns
    df = df[available_sort_columns]

    # Print a list of the best reference genomes in sorted order
    if sort_by in ["N50", "N90", "Length", "Score"]:
        df.sort_values(sort_by, inplace=True, ascending=False)
    else:
        df.sort_values(sort_by, inplace=True, ascending=True)
    df = df[0: top_n]
    df.to_csv(sys.stdout, index=False, sep="\t", float_format="%e")
