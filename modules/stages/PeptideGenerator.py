from __future__ import annotations

import logging
import re
import statistics
import sys
from dataclasses import dataclass, field
from io import StringIO
from operator import attrgetter
from typing import List

from VIPER import configmanager as cm

import pandas as pd

import custom_funcs


class REBprocessor:
    """
    Parses results from Residue Energy Breakdown
    """

    energy_cutoff: float = None
    max_side_extension: int = None
    max_length: int = None
    ignore_neighbors: bool = None
    always_include_direct_neighbors: bool = None
    prefilter_by_energy: bool = None
    max_neighbor_distance: int = None

    def __init__(self, energy_cutoff: float = -0.01, max_side_extension: int = 2, max_length: int = -1,
                 ignore_neighbors: bool = False, always_include_direct_neighbors: bool = False,
                 prefilter_by_energy: bool = False, max_neighbor_distance: int = 1):
        """
        Initializes the residue energy breakdown parser.

        :param energy_cutoff: Cutoff for which residue interactions to include, noninclusive. Any residue interactions
            with abs(total energy) < energy_cutoff will not be included. Default: 0.01
        """
        # TODO: Initialize through ConfigManager instead
        self.energy_cutoff = energy_cutoff
        self.max_side_extension = max_side_extension
        self.max_length = max_length
        self.ignore_neighbors = ignore_neighbors
        self.always_include_direct_neighbors = always_include_direct_neighbors
        self.prefilter_by_energy = prefilter_by_energy
        self.max_neighbor_distance = max_neighbor_distance

    def read_in(self, breakdown_file: str) -> tuple:
        """
        Parses the output from a residue energy breakdown run and prepares it for further analysis.

        :param breakdown_file:
        :return:
        """
        totals = []
        logging.info(f"Reading in residue energy breakdown file {breakdown_file}...")
        formatted = ""
        # Reformat to CSV-like
        with open(breakdown_file, "r") as i:
            pattern = re.compile(" +")
            for line in i:
                formatted += re.sub(pattern, ",", line)
        csv = pd.read_csv(StringIO(formatted),
                          usecols=lambda c: c.upper() in ["PDBID1", "RESTYPE1", "PDBID2", "RESTYPE2", "TOTAL"])

        # Get amino acids on chains and interactions per amino acid
        # chains = {}
        # interactions = {}
        nlist = []
        seen = {}
        for row in csv.itertuples(index=False, name="Interaction"):
            if (row.restype2 == "onebody") or (
                    self.prefilter_by_energy and (float(row.total) < self.energy_cutoff)):
                continue
            # Create Node and parse data
            if row.pdbid1 not in seen:
                node = REBprocessor.Node(amino_acid=row.restype1, residue=int(row.pdbid1[:-1]), chain=row.pdbid1[-1],
                                         partners=[(row.pdbid2, float(row.total))])
                if row.pdbid2[-1] in node.strength:
                    node.strength[row.pdbid2[-1]] += row.total
                else:
                    node.strength[row.pdbid2[-1]] = row.total
                seen[row.pdbid1] = node
                nlist.append(node)
                totals.append(row.total)
            else:  # Update Node with parsed data
                node = seen[row.pdbid1]
                node.partners.append((row.pdbid2, float(row.total)))
                totals.append(row.total)
                if row.pdbid2[-1] in node.strength:
                    node.strength[row.pdbid2[-1]] += row.total
                else:
                    node.strength[row.pdbid2[-1]] = row.total
            continue
        '''
            if row.restype2 == "onebody":
                if row.pdbid1[-1] in chains:
                    chains[row.pdbid1[-1]].append([row.pdbid1, row.restype1])
                else:
                    chains[row.pdbid1[-1]] = [[row.pdbid1, row.restype1]]
            else:
                totals.append(row.total)
                if row.pdbid1 in interactions:
                    interactions[row.pdbid1].append([row.pdbid1, row.pdbid2, row.total, row.restype1])
                    if row.pdbid2 in interactions:
                        interactions[row.pdbid2].append([row.pdbid1, row.pdbid2, row.total, row.restype1])
                    else:
                        interactions[row.pdbid2] = [[row.pdbid1, row.pdbid2, row.total, row.restype1]]
                else:
                    interactions[row.pdbid1] = [[row.pdbid1, row.pdbid2, row.total, row.restype1]]
                    if row.pdbid2 in interactions:
                        interactions[row.row.pdbid2].append([row.pdbid1, row.pdbid2, row.total, row.restype1])
                    else:
                        interactions[row.pdbid2] = [[row.pdbid1, row.pdbid2, row.total, row.restype1]]
        '''
        # Update neighbors
        for residueid, node in seen.items():
            prev_id = str(int(residueid[:-1]) - 1) + residueid[-1]
            if prev_id in seen:
                node.neighbor_prev = seen[prev_id]
            next_id = str(int(residueid[:-1]) + 1) + residueid[-1]
            if next_id in seen:
                node.neighbor_next = seen[next_id]

        stats = {"mean": statistics.mean(totals),
                 "stdev": statistics.stdev(totals),
                 "min": min(totals),
                 "max": max(totals)}
        return nlist, stats

    @dataclass
    class Node:
        _nodes = []
        amino_acid: str = None
        residue: int = None
        chain: str = None
        partners: List[tuple] = None  # (<partner residue id + chain id>, <interaction energy>)
        strength: dict = field(default_factory=lambda: ({}))
        neighbor_prev: REBprocessor.Node = None
        neighbor_next: REBprocessor.Node = None

        def __post_init__(self):
            REBprocessor.Node._nodes.append(self)

        # FIXME: This is a bad idea if multiple files are analyzed, can't be differentiated from
        #  where the multiple "19A" are from
        @classmethod
        def get(cls, residue_id: str) -> list:
            return [node for node in cls._nodes if residue_id in node]

        def __contains__(self, item):
            return item == str(self.residue) + self.chain

        def __repr__(self):
            return str(self.residue) + self.chain + " " + str(self.amino_acid)

    def transform_interactions(self, chain: str, interactions: dict) -> list:
        tally = []
        prev_store = None
        first = True
        for index, item in enumerate(interactions.items()):
            if item[0][-1] != chain:  # Ignore residues not on the submitted chain
                continue
            curr_node = None
            total_energy = 0  # This is only the total interaction energy with residues on _other_ chains!
            partners = []
            for interaction in item[1]:
                if interaction[1][-1].upper() == chain.upper():
                    continue
                total_energy += interaction[2]
                partners.append((interaction[1], interaction[2]))
            if len(partners) == 0:  # Ignore residues that only interact with their own chain
                continue
            if first:
                curr_node = self.Node(amino_acid=item[1][0][3], residue=item[0], partners=partners,
                                      strength=total_energy)
                first = False
            else:
                if abs(int(item[0][:-1]) - int(prev_store.residue[:-1])) <= self.max_neighbor_distance:
                    curr_node = self.Node(amino_acid=item[1][0][3], residue=item[0], partners=partners,
                                          strength=total_energy, neighbor_prev=prev_store)
                    tally[-1].neighbor_next = curr_node
                else:
                    curr_node = self.Node(amino_acid=item[1][0][3], residue=item[0], partners=partners,
                                          strength=total_energy)
            prev_store = curr_node
            tally.append(curr_node)
        return tally

    def reduce(self, from_chain: str, to_chain: str, nodes: list) -> list:

        def _include(add_to: list, n: REBprocessor.Node, curr_depth: int = 0) -> None:
            if False:  # cm.get("rosetta_config.custom_reb_func"):
                custom_funcs.rosetta_energy_breakdown_node_inclusion(add_to, n, depth=curr_depth, REBproc=self)
            else:
                if self.max_length == -1:
                    pass
                elif len(add_to) > self.max_length or to_chain not in n.strength:
                    return
                if self.always_include_direct_neighbors and curr_depth == 1:
                    add_to.append(n)
                if n.strength.get(to_chain, 100000) < self.energy_cutoff and curr_depth <= self.max_side_extension or (
                        self.always_include_direct_neighbors and curr_depth == 1):
                    if n not in add_to:
                        add_to.append(n)
                    if self.ignore_neighbors:
                        return
                    prev = None
                    nxt = None
                    if p := n.neighbor_prev:
                        prev = p
                    if nx := n.neighbor_next:
                        nxt = nx
                    if prev and nxt:
                        if n.neighbor_prev.strength.get(to_chain, 100000) < n.neighbor_next.strength.get(to_chain,
                                                                                                         100000):
                            _include(add_to, n.neighbor_prev, curr_depth + 1)
                            _include(add_to, n.neighbor_next, curr_depth + 1)
                        else:
                            _include(add_to, n.neighbor_next, curr_depth + 1)
                            _include(add_to, n.neighbor_prev, curr_depth + 1)
                    elif prev:
                        _include(add_to, n.neighbor_prev, curr_depth + 1)
                    elif nxt:
                        _include(add_to, n.neighbor_next, curr_depth + 1)

        final_residues = []
        elligible_nodes = [n for n in nlist if n.chain == from_chain and to_chain in n.strength]
        if len(elligible_nodes) == 0:
            logging.warning(
                f"Couldn't determine any residues that interacted significantly according to the specifications. Exiting...")
            sys.exit()
        for residue in sorted(elligible_nodes, key=lambda node: node.strength[to_chain], reverse=False):
            _include(final_residues, residue)
        return sorted(final_residues, key=lambda node: node.residue, reverse=False)


if __name__ == "__main__":
    p = REBprocessor()
    nlist, _ = p.read_in("energy_breakdown.out")
    # nlist = p.transform_interactions("A", nlist)

    result = p.reduce("A", "E", nlist)
    input()
