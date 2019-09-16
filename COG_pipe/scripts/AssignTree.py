#!/usr/bin/env python
import sys
import re

from ete3 import Tree
from collections import defaultdict, Counter
from operator import itemgetter

MIN_SUPPORT = 0.0


def read_lineage_file(lineage_file):
    mapping = {}
    for line in open(lineage_file):
        line = line.rstrip()
        (taxa_id, domain_id, phylum_id, class_id, order_id,
         family_id, genus_id, species_id) = line.split("\t")
        mapping[int(taxa_id)] = [domain_id, phylum_id, class_id,
                                order_id, family_id, genus_id, species_id]
    return mapping


def read_mapping_file(mapping_file):
    mapping = {}
    for line in open(mapping_file):
        line = line.rstrip()
        (taxa_id, name) = line.split(",")
        mapping[name] = int(taxa_id)
    return mapping


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("tree_file", help="Newick tree file")
    parser.add_argument("map_file", help="Maps names to lineages")
    parser.add_argument("lineage_file", help="Lineages")
    args = parser.parse_args()
    #import ipdb; ipdb.set_trace()
    t = Tree(args.tree_file)
    mapping = read_mapping_file(args.map_file)
    lineages = read_lineage_file(args.lineage_file)
    leaves = t.get_leaves()
    clusters = []
    references = []
    for leaf in leaves:
        if re.search(r".*Bin.*", leaf.name):
            clusters.append(leaf)
        else:
            references.append(leaf)
    # get closest matching reference for each cluster
    cluster_assign = defaultdict(dict)
    for cluster in clusters:
        closest, closest_dist = cluster.get_farthest_node()
    #    print cluster.name
     #   sys.stdout.flush()
        for reference in references:
            dist = t.get_distance(cluster, reference, topology_only=False)
            # print dist
            if dist < closest_dist:
                closest = reference
                closest_dist = dist
        ancestor = t.get_common_ancestor(cluster, closest)
        assigns = Counter()
        matches = ancestor.get_leaves()
        collate_hits = list()
        for depth in range(7):
            collate_hits.append(Counter())
        for match in matches:
            if match.name in mapping:
                taxa_id = mapping[match.name]
                if taxa_id in lineages:
                    line_match = lineages[taxa_id]
                    for depth in range(7):
                        collate_hits[depth][line_match[depth]] += 1.0
        for depth in range(7):
            collate = collate_hits[depth]
            d_weight = sum(collate.values())
            sort_collate = sorted(
                collate.items(), key=operator.itemgetter(1), reverse=True)
            nL = len(collate)
            if nL > 0:
                dp = 0.0
                if d_weight > 0.0:
                    dp = float(sort_collate[0][1])/d_weight
                    cluster_assign[cluster.name][depth] = (sort_collate[0][0], dp, d_weight)
            else:
                cluster_assign[cluster.name][depth] = ('None', 1.0, 1.0)
    for cluster in cluster_assign.keys():
        (assign, p, W) = cluster_assign[cluster][0]
        sys.stdout.write('%s' % cluster)
        if W >= MIN_SUPPORT:
            for depth in range(7):
                (assign, p, W) = cluster_assign[cluster][depth]
                sys.stdout.write('\t%d,%s,%f,%f' % (depth, assign, p, W))
            sys.stdout.write('\n')
            sys.stdout.flush()
        else:
            for depth in range(7):
                (assign, p, W) = cluster_assign[cluster][depth]
                sys.stdout.write('\t%d,NA,0.0,0.0' % (depth))
            sys.stdout.write('\n')
            sys.stdout.flush()


if __name__ == "__main__":
    main(sys.argv[1:])
