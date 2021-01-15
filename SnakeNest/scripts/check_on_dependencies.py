#!/usr/bin/env python
from subprocess import Popen, PIPE
from os.path import dirname,realpath,abspath
from collections import defaultdict
import sys 

def command_does_not_works(cmd,store_ans=False):
	command = Popen([cmd], stdout=PIPE, stderr=PIPE,shell=True)
	ans = command.communicate()[0]
	if store_ans:
		return command.returncode,ans
	else:
		return command.returncode


# this file should be placed at STRONG_dir/snakenest/scripts
STRONG_dir = dirname(dirname(dirname(abspath(realpath(__file__)))))

# check on R libraries
issues = defaultdict(list)
Rlib = [line.rstrip() for line in open("%s/R_libs.txt"%STRONG_dir)]
for lib in Rlib:
	if command_does_not_works("R -e 'library(%s)' "%lib):
		issues["Missing R library :"].append(lib)

# check on that annoying linkage for a R libraries
if command_does_not_works("ls $CONDA_PREFIX/lib/R/modules/libRlapack.so"):
	issues["R lapack symlink is not done :"].append("check the install readme for more info")

# check on desman
if command_does_not_works("desman -h"):
	issues["desman seems to be missing"].append("")

# check on spades
_,out = command_does_not_works("%s/SPAdes/assembler/build_spades/bin/unitig-coverage -h"%STRONG_dir,True)
if "GFA or prefix of the SPAdes binary saves" not in str(out):
	issues["SPAdes is not compiled"].append("")

# check on bayespath
if command_does_not_works("%s/BayesPaths/bin/bayespaths -h"%STRONG_dir):
	issues["bayespath seems to be missing"].append("")

# check on concoct refine
if command_does_not_works("concoct_refine -h"):
	issues["concoct_refine"].append("seems to be missing")

# check that concoct_refine has been fixed
_,found = command_does_not_works("grep 'int(NK), args.seed, args.threads)' $(which concoct_refine)",True)
if found:
	issues["concoct_refine"].append("need to be fixed : check the install readme for more info")



# coloration is done with https://askubuntu.com/questions/821157/print-a-256-color-test-pattern-in-the-terminal
if issues: 
	sys.exit("------------------------------------------------------------------------------\nFollowing issues where found :\n%s\n------------------------------------------------------------------------------"%"\n".join(["\033[38;5;196m%s %s\033[0m"%(issue_type,issue) for issue_type,issuelist in issues.items() for issue in issuelist]))
else:
	print('\033[0;32m no issue found with install\033[0m')












