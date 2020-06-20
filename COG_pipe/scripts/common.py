from __future__ import print_function
from __future__ import division

try:
    from future_builtins import zip
except:
    pass

import os
import os.path
import re

default_values = {
    "concoct_fragment_size": 10000,
    "concoct_contig_size": 500,
    "concoct_bin_multiplier": 3,
    "concoct_bin_max": 2000,
    "threads":     8,
    "assembly":    {"assembler": "spades", "k": [21, 33, 55],
                    "mem": 120, "threads": 16, "groups": []},
    "desman": {"execution": 0, "nb_haplotypes": 10, "nb_repeat": 5,
               "min_cov": 1, "dscripts": None},
    "bayespaths": {"min_orf_number_to_merge_bins": 10,
                   "min_orf_number_to_run_a_bin": 10,
                   "percent_unitigs_shared": 0.1,
                   "nb_strains": 16,
                   "nmf_runs":10,
                   "max_giter":4},
    "maganalysis": {"execution": 0},
    "evaluation": {"execution": 0, "genomes": ""},
}

# Taken from http://stackoverflow.com/questions/36831998/how-to-fill-default-parameters-in-yaml-file-using-python


def setdefault_recursively(tgt, default=default_values):
    for k in default:
        if isinstance(default[k], dict):  # if the current item is a dict,
            # expand it recursively
            setdefault_recursively(tgt.setdefault(k, {}), default[k])
        else:
            # ... otherwise simply set a default value if it's not set before
            tgt.setdefault(k, default[k])


def fill_default_values(config):
    local_dir = config.get("LOCAL_DIR")
    if local_dir:
        default_values["scripts"] = os.path.join(local_dir, "scripts")
        default_values["scg_data"] = os.path.join(local_dir, "scg_data")
        default_values["bayespaths"]["dir"] = os.path.join(
            local_dir, "..", "BayesPaths")
        default_values["desman"]["dscripts"] = os.path.join(
            local_dir, "..", "DESMAN/scripts")
        default_values["evaluation"]['scripts'] = os.path.join(
            local_dir, "scripts/evaluation")
        default_values["spades_tools"] = os.path.join(
            local_dir, "..", "SPAdes/assembler/build_spades/bin")
        default_values["assembly"]['dir'] = os.path.join(
            local_dir, "..", "SPAdes/assembler/bin")
    setdefault_recursively(config)


def sample_name(fullname):
    return os.path.splitext(os.path.basename(fullname))[0]


FASTA_EXTS = {".fasta", ".fasta.gz", ".fa", ".fna", ".fsa",
              ".fastq", ".fastq.gz", ".fq", ".fq.gz", ".fna.gz"}


def gather_paths(path, basename=False):
    for filename in os.listdir(path):
        name = os.path.basename(filename)
        for ext in FASTA_EXTS:
            if not name.endswith(ext):
                continue
            filepath = os.path.join(path, filename)
            if basename:
                yield (name[0:-len(ext)], filepath)
            else:
                yield filepath


def detect_reads(dir):
    return sorted(list(gather_paths(dir)))[:2]
