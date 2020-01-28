#!/usr/bin/env python3
import argparse
import subprocess
from subprocess import PIPE,Popen
import sys
import os
import os.path
import shutil
import yaml
import time

from scripts.common import fill_default_values

#copied from http://stackoverflow.com/questions/431684/how-do-i-cd-in-python/13197763#13197763
class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

parser = argparse.ArgumentParser(description="STRONG - STrain Resolution ON Graphs")
parser.add_argument("--threads", "-t", type=int, default=1, help="Number of threads")
parser.add_argument("dir", type=str, help="Output directory")
parser.add_argument("--config", "-c", type=str, default="", help="config.yaml to be copied to the directory (unnecessary if config.yaml is already there)")
parser.add_argument("--verbose", "-v", action="store_true", help="Increase verbosity level")
parser.add_argument("--dryrun", action="store_true", help="Show tasks, do not execute them")
parser.add_argument("--unlock", "-u", action="store_true", help="Unlock the directory")
parser.add_argument("--dag", "-d", help="file where you want the dag to be stored")
parser.add_argument('-s', nargs=argparse.REMAINDER,help="Pass additional argument directly to snakemake")
args = parser.parse_args()

exec_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
LOCAL_DIR = os.path.realpath(exec_dir)

base_params = ["snakemake", "--directory", os.path.realpath(args.dir), "--cores", str(args.threads), "--config", "LOCAL_DIR" + "=" + LOCAL_DIR, "--latency-wait", "120","-k"]

if args.verbose:
    # Output commands + give reasons + verbose (add "-n" for dry-run)
    base_params.extend(["-p", "-r", "--verbose"]) 
if args.dryrun:
    base_params.extend(["--dryrun"])
if args.unlock:
    base_params.extend(["--unlock"])
if not os.path.exists(args.dir):
    os.makedirs(args.dir)
if args.dag:
    base_params.extend(["--rulegraph"])
if args.s :
    base_params.extend(args.s)

print("Output folder set to", args.dir)

config_path = os.path.join(args.dir, "config.yaml")
if args.config:
    if os.path.exists(config_path):
        if subprocess.call(["diff", config_path, args.config]):
            print("Config path specified, but different config.yaml already exists in output folder", args.dir)
            sys.exit(239)
    else:
        print("Copying config from", args.config)
        shutil.copy(args.config, config_path)

with cd(exec_dir):
    call_cnt=0
    def call_snake(extra_params=[]):
        global call_cnt
        call_cnt+=1
        if args.dag:
            p1=Popen(base_params + extra_params, stdout=PIPE, stderr=sys.stderr)
            p2=Popen(["dot","-Tpng"], stdin=p1.stdout, stdout=PIPE, stderr=sys.stderr)
            with open(args.dag.replace(".png", str(call_cnt)+".png"), "bw") as f:
                f.write(p2.communicate()[0])
        else:
            subprocess.check_call(base_params + extra_params, stdout=sys.stdout, stderr=sys.stderr)

    def reuse_dir(dir_from, dir_name):
        if not dir_from:
            return
        local_dir = os.path.join(args.dir, dir_name)
        if not os.path.isdir(dir_from):
            print("Warning: {} source directory doesn't exist".format(dir_from))
            return
        if os.path.exists(local_dir):
            print("Warning: {} destination directory already exists".format(dir_name))
            return
        os.symlink(dir_from, local_dir)

    with open(config_path) as config_in:
        config = yaml.load(config_in)
    fill_default_values(config)
    
    print("Step #1 - Assembly / Binning / COG Annotation")
    call_snake(["--snakefile", "SCogSubGraph.snake"])


    print("Step #2 - Bin annotation")
    call_snake(["--snakefile", "Bin_annotation.snake"])

    #print("Step #2 - Subgraph Processing / Bin merging")
    #call_snake(["--snakefile", "HeavyLifting.snake"])

    #print("Step #3 - Strain Decomposition")
    #call_snake(["--snakefile", "BayesAGraph.snake"])

    if config["desman"]["execution"]:
        print("Step #4 - Strain Analysis with Desman") 
    #     #TODO (for Sergey) return the previously removed checkpoint and remove two calls 
    #     #call_snake(["--snakefile", "Desman.snake", "prepare"])
        #call_snake(["--snakefile", "Desman.snake", "all"])

    if config["maganalysis"]["execution"]:
        print("Step #5 - MAGAnalysis : place mags in a tree of references") 
        #call_snake(["--snakefile", "MAGAnalysis.snake"])
        
    if config["evaluation"]["execution"]:
        print("Step #6 - running evaluation") 
     #   call_snake(["--snakefile", "eval.snake"])
