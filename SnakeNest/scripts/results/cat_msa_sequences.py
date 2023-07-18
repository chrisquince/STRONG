#!/usr/bin/env python
# -*- coding: latin-1 -*-
import os
import glob
import argparse
import numpy as np 
from os.path import basename, dirname
from collections import Counter, defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp


is_mag = lambda name:(("_haplo_" not in name)&('Bin_' in name)) & (not "Hifi" in name)

def get_cog_coverage(scg_cov, map_name_to_orfs, cog_to_strain_seq):
    cog_mag_to_orf = {(line.rstrip().split("\t")[0].split(" ")[1],line.rstrip().split("\t")[1]):line.rstrip().split("\t")[0].split(" ")[0] for line in open(map_name_to_orfs)}
    orfs_to_covtot = {line.split("\t")[0]:sum([float(el) for el in line.split("\t")[1:]]) for line in open(scg_cov) if line.split("\t")[0] in set(cog_mag_to_orf.values())}
    cog_to_haplo_cov_tot = {cog:{haplo:orfs_to_covtot[cog_mag_to_orf[(cog,haplo)]] for haplo in haplo_seq if is_mag(haplo)} for cog,haplo_seq in cog_to_strain_seq.items()}
    return cog_to_haplo_cov_tot




def load_matrix(file,sample_order=None,strain_order=None) :
    with open(file) as handle :
        header = next(handle).rstrip().split("\t")[1:]
        strains = []
        matrix = []
        for line in handle : 
            # if there are non overlapping sequences it can yield NA
            splitlines = line.rstrip().replace("NA","0").split("\t")
            strains.append(splitlines[0])
            matrix.append(list(map(float,splitlines[1:])))
    matrix = np.array(matrix)
    if sample_order :
        reorder_samples = [header.index(sample) for sample in sample_order]
        reorder_strain = [strains.index(strain) for strain in strain_order]
        return matrix[:,reorder_samples][reorder_strain,:]
    else : 
        return matrix,header,strains


def similarity(mat):
    # does not make much sense, percent substitution is also a good mesure of distance,
    # but precision is not good, looking at similarity is a dumb way to tackle ambiguity (multiple min dist)
    s = np.zeros(mat.shape)
    for row,i in enumerate(mat):
        for col,j in enumerate(mat):
            Ni = np.linalg.norm(i)
            Nj = np.linalg.norm(j)
            if Ni==0:
                i+=1e-10
                Ni = np.linalg.norm(i)
            if Nj==0:
                j+=1e-10
                Nj = np.linalg.norm(j)
            iN = i/np.linalg.norm(i)
            jN = j/np.linalg.norm(j)
            s[row,col] = np.dot(iN,jN)
    return s

def delete_seq(dist_mat,index_mag,rownames,colnames):
    # need to delete 1 by 1 
    identicals = (dist_mat[:,index_mag]==0).sum(1) 
    while max(identicals[index_mag])!=1:
        index_to_del = next(index for index,val in enumerate(identicals) if (index_mag[index])&(val==(max(identicals[index_mag]))))
        to_del = np.array([index==index_to_del for index in range(len(identicals))])
        reselect = ~(to_del*index_mag)
        dist_mat = dist_mat[reselect,:][:,reselect]
        colnames = rownames = list(np.array(rownames)[reselect])
        index_mag = index_mag[reselect]
        identicals = (dist_mat[:,index_mag]==0).sum(1)
    return dist_mat,index_mag,rownames,colnames


def remove_selected(haplo_cog_mags,cog_to_dists,selected):
    haplo_cog_mags = {hap:{cog:[name for name in names if name not in selected[cog]] for cog,names in cog_to_names.items()} for hap,cog_to_names in haplo_cog_mags.items()}
    haplo_cog_mags = {hap:{cog:names for cog,names in cog_to_names.items() if names} for hap,cog_to_names in haplo_cog_mags.items()}
    haplo_cog_mags = {hap:cog_to_names for hap,cog_to_names in haplo_cog_mags.items() if cog_to_names!={}}
    cog_to_dists = {cog:{name:dist for name,dist in name_to_dist.items() if name not in selected[cog]} for cog,name_to_dist in cog_to_dists.items()}
    cog_to_dists = {cog:name_to_dist for cog,name_to_dist in cog_to_dists.items() if name_to_dist}
    return haplo_cog_mags, cog_to_dists

def select_best_seq(mag, hapo_ref,sorted_haplo, cog_to_dists, new_mag_to_old, cog_to_haplo_cov_tot = None):
    selected = defaultdict(list)
    for cog,name_dists in cog_to_dists.items():
        if hapo_ref == "cov":
            # filter par name_dists, because I can't be arsed to remove non pertinant seqs from cog_to_haplo_cov_tot like I did to dist_mat with delete_seq()
            candidate_main = [[name,cov] for name,cov in cog_to_haplo_cov_tot[cog].items() if name in name_dists]
        else : 
            candidate_main = [[name,dists[sorted_haplo.index(hapo_ref)]] for name,dists in name_dists.items()]
        select = max(candidate_main,key=lambda x:x[1])[0]
        new_mag_to_old[mag][cog] = select
        selected[cog].append(select)
    return selected

def get_new_name(mag):
    if "_nb_" in mag :
        name,counter = mag.split("_nb_")
        counter = int(counter)+1
        mag = "%s_nb_%s"%(name,counter)
    else : 
        mag = mag+"_nb_2"
    return mag



def mags_seqs(path,sorted_cogs,sorted_strains,cog_to_haplo_cov_tot):
    mags_names = [name for name in sorted_strains if is_mag(name)]
    sorted_haplo = sorted([name for name in sorted_strains if not is_mag(name)])
    if len(mags_names)>1:
        haplo_cog_mags = defaultdict(lambda:defaultdict(list))
        cog_haplo_mags = defaultdict(lambda:defaultdict(list))
        cog_to_dists = defaultdict(lambda:{})
        for cog in sorted_cogs:
            file = "%s/%s_dist_mat.tsv"%(path,cog)
            if open(file).read()=="":
                # sometimes there is only 1 seq for a cog, we could use a checkpoint, to just no do this, but it's super heavy dag wise. 
                # so let's just no do that on go on when a file is empty
                continue
            dist_mat,colnames,rownames = load_matrix(file)
            index_mag = np.array([is_mag(name) for name in rownames])

            # in some case dist is nan from fragmentation, let's just say it's dist 0 and remove them as identicals
            dist_mat[np.isnan(dist_mat)] = 0

            # ignore any sequence which are identical to others, it can be due to fragmentation, or just the same seq multiple time
            dist_mat,index_mag,rownames,colnames = delete_seq(dist_mat,index_mag,rownames,colnames)

            # get similarity from distance, not super useful
            simil_mat = similarity(dist_mat)
            np.fill_diagonal(simil_mat,0)
            # get for each mag the best haplo it looks like
            for index,line in enumerate(simil_mat):
                if index_mag[index]:
                    haplo_best_fit = np.array(colnames)[(~index_mag)*(line==max(line[~index_mag]))][0] 
                    haplo_cog_mags[haplo_best_fit][cog].append(colnames[index])
                    cog_haplo_mags[cog][haplo_best_fit].append([colnames[index],max(line)])
                    cog_to_dists[cog][colnames[index]] = [line[colnames.index(strain)] if  strain in colnames else 0 for strain in sorted_haplo]

        # get for the main mag, choose the one with most coverage 
        new_mag_to_old = defaultdict(lambda:defaultdict(list))
        mag = mags_names[0].split("_nb_")[0]
        selected = select_best_seq(mag, "cov",sorted_haplo, cog_to_dists, new_mag_to_old, cog_to_haplo_cov_tot)
        haplo_cog_mags,cog_to_dists = remove_selected(haplo_cog_mags, cog_to_dists, selected)
        # all others, just create 1 seq by haplotype left
        while haplo_cog_mags!={}:
            for haplo in haplo_cog_mags:
                mag = get_new_name(mag)
                selected = select_best_seq(mag, haplo,sorted_haplo, cog_to_dists, new_mag_to_old)
                haplo_cog_mags,cog_to_dists = remove_selected(haplo_cog_mags, cog_to_dists, selected)
    else :
        mag = mags_names[0]
        new_mag_to_old = {mag:{cog:mag for cog in sorted_cogs}}
    return new_mag_to_old





def main(scg_cov, map_name_to_orfs, output):
    path = dirname(output)
    msa_files = glob.glob("%s/*_trim.msa"%path)
    dist_mat = glob.glob("%s/*_dist_mat.tsv"%path)

    strain_to_cog_seq = defaultdict(lambda:{})
    cog_to_strain_seq = defaultdict(lambda:{})
    for file in msa_files:
        cog = basename(file).replace("_trim.msa","")
        for strain,seq in sfp(open(file)):
            strain_to_cog_seq[strain][cog] = seq
            cog_to_strain_seq[cog][strain] = seq
    sorted_cogs = sorted(cog_to_strain_seq.keys())
    sorted_strains = sorted(strain_to_cog_seq.keys())
    mean_cog_len = {cog:int(np.mean([len(seq) for seq in dict_strain.values()])) for cog,dict_strain in cog_to_strain_seq.items()}

    ### deal with multiples mag cog : find which should go with which
    # get coverage tot of each mag cog orf
    cog_to_haplo_cov_tot = get_cog_coverage(scg_cov, map_name_to_orfs, cog_to_strain_seq)
    # choose the consensus sequence orfs, this way too complicated for what it does
    new_mag_to_old = mags_seqs(path,sorted_cogs,sorted_strains,cog_to_haplo_cov_tot)

    # get strain concatenated sequence
    sorted_strains = sorted([strain for strain in sorted_strains if not is_mag(strain)]+list(new_mag_to_old.keys()))
    strain_to_seq = defaultdict(str)
    for strain in sorted_strains:
        for cog in sorted_cogs: 
            if is_mag(strain):
                if cog in new_mag_to_old[strain]:
                    cog_to_seq = strain_to_cog_seq[new_mag_to_old[strain][cog]]
                else:
                    cog_to_seq = {}
            else:
                cog_to_seq = strain_to_cog_seq[strain]
            if cog in cog_to_seq:
                strain_to_seq[strain]+=cog_to_seq[cog]
            else :
                # deal with missing cogs
                strain_to_seq[strain]+=mean_cog_len[cog]*"-"
    with open(output,"w") as handle : 
        handle.writelines(">%s\n%s\n"%(strain,seq) for strain,seq in strain_to_seq.items())



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('scg_cov', help=("path to scg orfs coverage files"))
    parser.add_argument('map_name_to_orfs', help=("path to file mapping mag sequence names to orfs"))
    parser.add_argument('output', help=("output directory"))
    args = parser.parse_args()
    main(args.scg_cov, args.map_name_to_orfs, args.output)

