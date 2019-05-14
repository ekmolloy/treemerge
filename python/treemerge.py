
""""
This file is a python prototype of TreeMerge from the paper:

Molloy, E.K., Warnow, T. (2019). TreeMerge: A new method for improving the
scalability of species tree estimation methods.

Copyright (c) 2018-19 Erin K. Molloy
All rights reserved.

License: 3-Clause BSD,
see https://opensource.org/licenses/BSD-3-Clause
"""

import os
import os.path
import sys

scriptpath = os.path.realpath(__file__)
sys.path.append(scriptpath.rsplit("/", 1)[0])

import njmergepair
import njmerge2

import argparse
from copy import deepcopy
import dendropy
import networkx
import numpy
import subprocess
import time


def get_base_name(file_name, remove_extension=True, remove_path=True):
    x = file_name
    if remove_path:
        x = x.rsplit('/', 1)  # Remove path
        if len(x) > 1:
            x = x[1]
        else:
            x = x[0]
    if remove_extension:
        x = x.rsplit('.', 1)[0]  # Remove file extension
    return x


def name_treepair_file(path, filei, filej):
    basei = get_base_name(filei)
    basej = get_base_name(filej)
    x = "njmergepair-" + basei + "-and-" + basej + ".tre"
    if path != "":
        x = path + "/" + x
    return x


def name_nexspair_file(path, filei, filej):
    basei = get_base_name(filei)
    basej = get_base_name(filej)
    x = "nexuspair-" + basei + "-and-" + basej + ".txt"
    if path != "":
        x = path + "/" + x
    return x


def add_branch_lengths_with_paup(paup, pdm, tree, nexfile, outfile):
    """Add non-negative branch lengths to tree using PAUP*

    Parameters
    ----------
    paup : str
        PAUP* binary name
    pdm : dendropy phylogenetic distance matrix object
    tree : dendropy tree object
    nexfile : str
        output nexus file name
    outfile : str
        output tree file name

    Returns
    -------
    name_map : dict

    """
    taxa = [x.label for x in pdm._taxon_phylogenetic_distances]
    ntax = len(taxa)

    for l in tree.leaf_nodes():
        x = l.taxon.label
        l.taxon.label = "paupsafe" + x

    with open(nexfile, 'w') as f:
        f.write("#NEXUS\n")
        f.write("begin taxa;\n")
        f.write("    dimensions ntax=%d;\n" % ntax)
        f.write("    taxlabels ")
        f.write("paupsafe" + " paupsafe".join(taxa))
        f.write(";\n")
        f.write("end;\n")
        f.write("begin distances;\n")
        f.write("     format triangle=both;\n")
        f.write("     matrix\n")
        for xi in pdm.taxon_namespace:
            row = []
            for xj in pdm.taxon_namespace:
                if xi == xj:
                    row.append("0.000000")
                else:
                    dij = pdm._taxon_phylogenetic_distances[xi][xj]
                    row.append(str("%1.6f" % dij))
            f.write("          " + "paupsafe" + xi.label + " " +
                    " ".join(row) + "\n")
        f.write("          ;\n")
        f.write("end;\n")
        f.write("begin paup;\n")
        f.write("     set criterion=distance;\n")
        f.write("     dset negbrlen=prohibit objective=lsfit;\n")
        f.write("     constraints fulltree=")
        f.write(tree.as_string(schema="newick"))
        f.write("     nj constraints=fulltree enforce=yes brlens=yes;\n")
        f.write("     savetrees file=" + outfile + " brlens replace;\n")
        f.write("     q;\n")
        f.write("end;\n")

    cmd = paup + " -n " + nexfile

    try:
        with open(os.devnull, "w") as out:
            p = subprocess.call(cmd.split(), stdout=out, stderr=out)
            if p:
                raise Exception("System call " + cmd +
                                " had non-zero exit status\n")
    except (OSError, IOError) as e:
        raise Exception("System call " + cmd +
                        " was unable to be performed!\n")

    subprocess.call(["rm", nexfile])


def set_zeros_to_half_shortest(tree):
    """Set branches with length zero to 1/2 shortest branch

    Parameters
    ----------
    tree : dendropy tree object

    Returns
    -------
    tree : dendropy tree object

    """
    elens = [e.length for e in tree.postorder_edge_iter()]
    elens = filter(None, elens)
    elens = numpy.array(elens)

    min_elen = numpy.min(elens)
    half_min_elen = min_elen / 0.5

    for e in tree.postorder_edge_iter():
        if e.length == 0.0:
            e.length = half_min_elen

    return tree


def get_backbone_tree(tree1, tree2):
    """Constain tree1 to its shared leaf set with tree2

    Parameters
    ----------
    tree1 : dendropy tree object
    tree2 : dendropy tree object

    Returns
    -------
    tree1 : dendropy tree object

    """
    tree1 = deepcopy(tree1)

    leaves1 = njmergepair.get_leaf_set(tree1)
    leaves2 = njmergepair.get_leaf_set(tree2)
    shared = list(leaves1.intersection(leaves2))

    taxa = dendropy.TaxonNamespace(shared)

    tree1.retain_taxa_with_labels(shared)
    tree1.migrate_taxon_namespace(taxa)
    tree1.encode_bipartitions()

    return tree1


def map_splits_to_node_list(big_tree, lil_tree):
    """Map a split in little tree to a list of nodes in big tree

    NOTE: Because little tree is *contained* within big tree, more
          than one node in big tree can be mapped to the same split.

    Parameters
    ----------
    big_tree : dendropy tree object
    lil_tree : dendropy tree object

    Returns
    -------
    split_to_node_list : python dictionary
        keys - splits encoded as integers (read below!)
        values - lists of nodes in a dendropy tree object

    """
    lil_leafset = njmergepair.get_leaf_set(lil_tree)

    split_to_node_list = {}
    for node in big_tree.postorder_node_iter():
        big_leafset = njmergepair.get_leaf_set(node)
        shared_leafset = list(big_leafset.intersection(lil_leafset))

        if len(shared_leafset) > 0:
            split = lil_tree.taxon_namespace.taxa_bitmask(labels=shared_leafset)
            try:
                split_to_node_list[split].append(node)
            except KeyError:
                split_to_node_list[split] = []
                split_to_node_list[split].append(node)

    return split_to_node_list


def combine_two_trees_via_dscm(tree_AB, tree_BC):
    """Combines two trees via distance-based strict consensus merger

    A is the subset of leaves only in tree AB
    C is the subset of leaves only in tree BC
    B is the subset of leaves in both tree AB and tree BC

    Tree AB and tree BC must be equivalent on their shared leaf set B

    Parameters
    ----------
    tree_AB : dendropy tree object
    tree_BC : dendropy tree object

    Returns
    -------
    tree_BC : dendropy tree object
        tree BC with leaves from set A added

    """
    data = list(njmergepair.get_leaf_set(tree_AB).intersection(njmergepair.get_leaf_set(tree_BC)))
    taxa = dendropy.TaxonNamespace(data)

    # [incompatible] = are_two_trees_incompatible(tree_AB, tree_BC)
    # if incompatible:
    #     sys.exit("Input trees are not compatible!\n")

    backbone_tree = get_backbone_tree(tree_AB, tree_BC)
    if backbone_tree is None:
        raise Exception("Unable to extract a backbone tree!\n")

    # Root all trees at the same shared leaf --
    # required for split mapping and post-order traversal to work!
    root = backbone_tree.taxon_namespace[0].label

    node_XX = backbone_tree.find_node_with_taxon_label(root)
    node_AB = tree_AB.find_node_with_taxon_label(root)
    node_BC = tree_BC.find_node_with_taxon_label(root)

    elen_AB = node_AB.edge.length / 2.0
    elen_BC = node_BC.edge.length / 2.0

    backbone_tree.is_rooted = True
    tree_AB.is_rooted = True
    tree_BC.is_rooted = True

    backbone_tree.reroot_at_edge(node_XX.edge)
    tree_AB.reroot_at_edge(node_AB.edge)
    tree_BC.reroot_at_edge(node_BC.edge)

    node_AB.edge.length = elen_AB
    node_BC.edge.length = elen_BC

    node_AB.sibling_nodes()[0].edge.length = elen_AB
    node_BC.sibling_nodes()[0].edge.length = elen_BC

    # Map nodes based on splits in shared leaf set
    map_AB = map_splits_to_node_list(tree_AB, backbone_tree)
    map_BC = map_splits_to_node_list(tree_BC, backbone_tree)

    # Add missing taxa from AB **into** BC using the
    # distance-based SCM strategy to handle collisions
    nodes = [n for n in backbone_tree.postorder_node_iter()]
    for node in nodes[:-1]:
        clade = njmergepair.get_leaf_set(node)
        split = backbone_tree.taxon_namespace.taxa_bitmask(labels=clade)

        node_path_AB = map_AB[split]
        node_path_BC = map_BC[split]

        num_edges_AB = len(node_path_AB)
        num_edges_BC = len(node_path_BC)

        sibs_path_AB = []
        for n in node_path_AB:
            s = n.sibling_nodes()
            if len(s) > 1:
                sys.exit("Tree AB is not binary!\n")
            sibs_path_AB.append(s[0])

        sibs_path_BC = []
        for n in node_path_BC:
            s = n.sibling_nodes()
            if len(s) > 1:
                sys.exit("Tree BC is not binary!\n")
            sibs_path_BC.append(s[0])

        if num_edges_AB > 1 and num_edges_BC > 1:
            # Found a collision -- add edges from tree AB to tree BC
            # Find normalization factor for path in AB
            # ALSO compute the point of attachment for edges in the path
            # IN ORDER to identify the ORDER in which edges should be added
            elen_AB = node_path_AB[0].edge.length
            elen_path_AB = [elen_AB]
            i_AB = 1
            for node_AB in node_path_AB[1:]:
                elen_AB = elen_AB + node_AB.edge.length
                elen_path_AB.append(elen_path_AB[i_AB - 1]
                                    + node_AB.edge.length)
                i_AB = i_AB + 1

            elen_BC = node_path_BC[0].edge.length
            elen_path_BC = [elen_BC]
            i_BC = 1
            for node_BC in node_path_BC[1:]:
                elen_BC = elen_BC + node_BC.edge.length
                elen_path_BC.append(elen_path_BC[i_BC - 1]
                                    + node_BC.edge.length)
                i_BC = i_BC + 1

            if elen_AB == 0.0 or elen_BC == 0.0:
                raise Exception("Collision on path of length zero!\n")

            norm_AB = elen_BC / elen_AB

            for i in range(len(elen_path_AB)):
                elen_path_AB[i] = elen_path_AB[i] * norm_AB

            # Extract components of tree BC
            node_BC = node_path_BC[-1]
            sibs_BC = sibs_path_BC[-1]
            parent_BC = node_BC.parent_node
            parent_BC.clear_child_nodes()

            # Get node in tree BC and update branch length
            child1 = node_path_BC[0]
            elen_AB = node_path_AB[0].edge.length * norm_AB
            elen_BC = node_path_BC[0].edge.length
            if elen_AB < elen_BC:
                child1.edge.length = elen_AB
            else:
                child1.edge.length = elen_BC

            i_AB = 0
            i_BC = 0
            start = None
            while (True):
                dothis = None
                if i_AB < num_edges_AB - 1 and i_BC < num_edges_BC - 1:
                    if elen_path_AB[i_AB] == elen_path_BC[i_BC]:
                        # Add a new node created from AB and BC
                        dothis = 3
                    elif elen_path_AB[i_AB] < elen_path_BC[i_BC]:
                        # Add remaining edges from tree AB
                        dothis = 1
                    else:
                        # Add remaining edges from tree BC
                        dothis = 2
                elif i_AB < num_edges_AB - 1:
                    # Add remaining edges from tree AB
                    dothis = 1
                elif i_BC < num_edges_BC - 1:
                    # Add remaining edges from tree BC
                    dothis = 2
                else:
                    # No more edges to add!
                    break

                if dothis == 1:
                    # Adding AB only
                    child2 = sibs_path_AB[i_AB]
                    if start is None:
                        start = child1.edge.length
                    else:
                        stop = elen_path_AB[i_AB]
                        child1.edge.length = stop - start
                        start = stop
                    i_AB = i_AB + 1
                elif dothis == 2:
                    # Add BC only
                    child2 = sibs_path_BC[i_BC]
                    if start is None:
                        start = child1.edge.length
                    else:
                        stop = elen_path_BC[i_BC]
                        child1.edge.length = stop - start
                        start = stop
                    i_BC = i_BC + 1
                else:
                    # Add both AB and BC
                    child2 = dendropy.Node()
                    child2.set_child_nodes([sibs_path_AB[i_AB],
                                            sibs_path_BC[i_BC]])
                    child2.edge.length = 0.0
                    if start is None:
                        start = child1.edge.length
                    else:
                        stop = elen_path_BC[i_BC]
                        child1.edge.length = stop - start
                        start = stop
                    i_AB = i_AB + 1
                    i_BC = i_BC + 1

                child2.parent_node = None

                # Combine nodes from tree BC (node 1) and
                # tree AB or tree BC (node 2)
                next_node = dendropy.Node()
                next_node.set_child_nodes([child1, child2])

                # Set node in tree BC as next_node
                child1 = next_node

            # Recombine the three components of tree BC
            next_node.edge.length = elen_path_BC[-1] - start
            parent_BC.set_child_nodes([next_node] + [sibs_BC])
        elif num_edges_AB > 1:
            # Found edges in tree AB not in tree BC --
            # add edges from tree AB to tree BC!

            # Find normalization factor for path in AB
            elen_BC = node_path_BC[0].edge.length
            elen_AB = 0.0
            elen_path_AB = []
            for node_AB in node_path_AB:
                elen_AB = elen_AB + node_AB.edge.length
                elen_path_AB.append(node_AB.edge.length)

            if elen_AB == 0.0:
                xxxx_AB = elen_BC / num_edges_AB
                for i in range(num_edges_AB):
                    elen_path_AB[i] = xxxx_AB
            else:
                norm_AB = elen_BC / elen_AB
                for i in range(num_edges_AB):
                    elen_path_AB[i] = elen_path_AB[i] * norm_AB

            # Extract components of tree BC
            node_BC = node_path_BC[0]
            sibs_BC = sibs_path_BC[0]
            parent_BC = node_BC.parent_node
            parent_BC.clear_child_nodes()

            # Get node in tree BC and update branch length to match tree AB
            child1 = node_path_BC[0]
            child1.edge.length = elen_path_AB[0]

            # Add each nodes in tree AB to tree BC
            for i_AB in range(1, num_edges_AB):
                # Remove parent from child 1
                child1.parent_node = None

                # Get node being added from tree AB!
                child2 = sibs_path_AB[i_AB-1]
                child2.parent_node = None

                # Combine nodes from tree BC (node 1) and tree AB (node 2)
                next_node = dendropy.Node()
                next_node.set_child_nodes([child1, child2])
                next_node.edge.length = elen_path_AB[i_AB]

                # Set node in tree BC as next_node
                child1 = next_node

            # Recombine the three components of tree BC
            parent_BC.set_child_nodes([next_node] + [sibs_BC])
        elif num_edges_BC > 1:
            # Found edges in tree BC not in tree AB
            pass
        else:
            # Found one edge in tree BC and one edge in tree AB
            pass

    tree_BC.migrate_taxon_namespace(taxa)


def dscmcombine(workdir, trees, mstmat, outfile):
    """
    Parameters
    ----------
    workdir : str
        working or output directory

    Results
    -------
    Nothing

    """
    # Turn matrix into graph
    graph = networkx.Graph(mstmat)

    # Merge trees
    combined_tree = None
    nodes = list(deepcopy(graph.nodes()))
    root = None
    next_roots = [nodes[0]]
    while len(nodes) > 0:
        if root is None:
            root = next_roots[0]
        neighbors = graph.neighbors(root)
        next_roots = list(set(next_roots).union(set(neighbors)))
        if len(neighbors) == 0:
            nodes.remove(root)
            next_roots.remove(root)
            root = None
        else:
            sys.stdout.write("Combining %d and %d...\n"
                             % (root, neighbors[0]))

            if root < neighbors[0]:
                i = root
                j = neighbors[0]
            else:
                i = neighbors[0]
                j = root

            tijfile = name_treepair_file(workdir, trees[i], trees[j])
            if not os.path.exists(tijfile):
                tijfile = name_treepair_file(workdir, trees[j], trees[i])

            if combined_tree is None:
                combined_tree = dendropy.Tree.get(path=tijfile,
                                                  schema="newick")
            else:
                tij = dendropy.Tree.get(path=tijfile, schema="newick")
                combine_two_trees_via_dscm(tij, combined_tree)
                combined_tree.update_bipartitions()

            sys.stdout.write("...combined tree has %d leaves!\n"
                             % len(njmergepair.get_leaf_set(combined_tree)))

            graph.remove_edge(root, neighbors[0])

    with open(outfile, 'w') as f:
        f.write(combined_tree.as_string(schema="newick")[5:])



def main(args):
    output = args.output
    if os.path.isfile(output):
        sys.exit(output + " already exists!\n")

    sys.stdout.write("Running TreeMerge...\n")
    total_ts = time.time()

    # Input
    paup = args.paup
    strefile = args.start
    treefiles = args.trees
    dmatfile = args.matrix
    taxafile = args.taxa
    workdir = args.workdir

    # Build MST from starting tree
    ts = time.time()
    mst = njmerge2.tree_to_mst(strefile, treefiles)
    rt = time.time() - ts
    sys.stdout.write("...computed MST in %d seconds.\n" % rt)

    # Merge pairs of trees uing NJMerge
    graph = networkx.Graph(mst)
    trees = []
    for e in graph.edges():
        ts = time.time()
        i, j = e

        tifile = treefiles[i]
        tjfile = treefiles[j]
        nexfile = name_nexspair_file(workdir, tifile, tjfile)
        tijname = name_treepair_file("", tifile, tjfile)
        tijfile = workdir + "/" + tijname

        ti = dendropy.Tree.get(path=tifile, schema="newick")
        tj = dendropy.Tree.get(path=tjfile, schema="newick")

        # Merge two trees using NJMerge
        [dij, tij] = njmergepair.run(dmatfile, taxafile, ti, tj)

        # Add (non-negative) branch lengths with PAUP*
        tij.is_rooted = False
        tij.collapse_basal_bifurcation(set_as_unrooted_tree=True)
        add_branch_lengths_with_paup(paup, dij, tij, nexfile, tijname)
        try:
            with open(tijfile, 'r'):
                pass
        except FileNotFoundError as e:
            raise Exception(tijfile + " was not created - " +
                            "check if PAUP* binary has expired!\n")

        # Set branches with length 0 to 1/2 shortest branch
        tij = dendropy.Tree.get(path=tijfile, schema="nexus")
        tij = set_zeros_to_half_shortest(tij)

        for l in tij.leaf_nodes():
            l.taxon.label = l.taxon.label.replace("paupsafe", "")

        ## Write tree
        with open(tijfile, 'w') as f:
            f.write(tij.as_string(schema="newick"))

        rt = time.time() - ts
        sys.stdout.write("...merged trees %d and %d in %d seconds.\n" % (i, j, rt))

    # Combine using the pairwise merged trees using distance-SCM to handle collisions
    ts = time.time()
    final_tree = dscmcombine(workdir, treefiles, mst, output)
    rt = time.time() - ts
    sys.stdout.write("...combined all trees in %d seconds.\n" % rt)

    total = time.time() - total_ts
    sys.stdout.write("Total execution time: %d seconds\n" % total)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--start", type=str,
                        help="Input starting tree file", required=True)

    parser.add_argument("-t", "--trees", type=str, nargs="+",
                        help="Input subset tree files", required=True)

    parser.add_argument("-m", "--matrix", type=str,
                        help="Input distance matrix", required=True)

    parser.add_argument("-x", "--taxa", type=str,
                        help="Taxon names corresponding to rows in distance matrix, necessary when using output from FastME",
                        required=False)

    parser.add_argument("-o", "--output", type=str,
                        help="Output tree",
                        required=True)

    parser.add_argument("-p", "--paup", type=str,
                        help="Path to PAUP* binary",
                        required=True)

    parser.add_argument("-w", "--workdir", type=str,
                        help="Path to working directory for writing temporary files",
                        required=False)

    main(parser.parse_args())
