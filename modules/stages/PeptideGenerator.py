from __future__ import annotations

import logging
import pprint
import sys
from abc import abstractmethod, ABCMeta
from dataclasses import dataclass
from math import sqrt
from typing import List, final, Type, Tuple

import ConfigManager
from modules.wrappers.RosettaWrapper import REBprocessor
from util import PDBtool

cm = ConfigManager.ConfigManager.get_instance


class PeptideGenerator:
    # https://new.rosettacommons.org/docs/latest/application_documentation/design/pepspec ?
    # Do multiple runs of relax + minimize
    #  Which residues consistently interact the strongest?
    # Use constraints file to pin those residues (movemap?)
    # Resfile PIKAA + amino acid letter
    #
    # https://www.rosettacommons.org/docs/latest/application_documentation/design/mpi-msd
    pass


class _SelectionStrategies:

    @staticmethod
    def add_linkers(nlist: List[REBprocessor.Node], respect_length_limit: bool = False) -> List[REBprocessor.Node]:
        """
        Adds linker nodes to every gap in consecutive residue ids in a node list and potentially also updates residue
        numbering, if necessary.

        :param nlist: The node list (candidate peptide) to add linkers to
        :param respect_length_limit: (Optional) Whether to raise an error if adding linkers would make the peptide
            length exceed the limit. (Default: False)
        :return: The original node list, but with linker nodes added and potentially updated numbering
        """
        linker = cm().get("peptide_generator.linker")
        if linker is None:
            linker = "GG"  # default is polyglycine (2)
        first = True
        new_nlist = []
        increase_orig_num = 0
        for n, node in enumerate(nlist):
            if first:
                new_nlist.append(REBprocessor.Node(amino_acid=node.amino_acid,
                                                   residue_id=node.residue_id,
                                                   chain=node.chain,
                                                   partners=node.partners,
                                                   strength=node.strength,
                                                   orig_res_id=node.orig_res_id))
                first = False
                continue
            interres = abs(node.residue_id - new_nlist[-1].residue_id)
            if interres > 1:  # Jump in the sequence, need linker
                # Reduce the original residue offset by the length of the gap, since all following original residues
                # will have their id increased by at least the length of the gap
                increase_orig_num = max(0, increase_orig_num - interres)
                # Number of elements to use from linker (usually it's all of them)
                use_num = len(linker)
                if interres - 1 < len(linker):
                    logging.warning(f"Trying to link a gap between fragments ({new_nlist[-1]} -> {node}) that is "
                                    f"shorter than the linker ({linker})! "
                                    f"Using policy {cm().get('peptide_generator.linker_oversize_policy').upper()}.")
                    policy = cm().get('peptide_generator.linker_oversize_policy').upper()
                    if policy == "TRUNCATE":
                        use_num = interres - 1
                    elif policy == "IGNORE":
                        # If we're inserting the full length linker anyway, we have to increase the following residue
                        # ids accordingly
                        increase_orig_num += use_num
                    elif policy == "SKIP":
                        # Don't add linker
                        # Copy over original node with (potentially) adjusted residue id
                        add_node = REBprocessor.Node(amino_acid=node.amino_acid,
                                                     residue_id=node.residue_id + increase_orig_num,
                                                     chain=node.chain,
                                                     partners=node.partners,
                                                     strength=node.strength,
                                                     neighbor_prev=new_nlist[-1],
                                                     orig_res_id=node.orig_res_id)
                        new_nlist[-1].neighbor_next = add_node
                        new_nlist.append(add_node)
                        continue
                    else:  # Default case: TRUNCATE
                        use_num = interres - 1
                base_id = new_nlist[-1].residue_id

                # Check that length of current candidate w/ linkers + linker + remaining peptide doesn't exceed length
                # (Doesn't check for future linkers, but they will be considered once they are encountered)
                if respect_length_limit and len(new_nlist) + len(linker) + len(nlist[n:]) >= cm().get(
                        "peptide_generator.max_length"):
                    logging.warning(f"Cannot insert linker, because it would make the peptide length exceed the limit. "
                                    f"Returning early with {new_nlist}")
                    return new_nlist

                for lnum, element in enumerate(linker, start=1):  # Add linkers
                    if lnum > use_num:  # stop early, because of truncate limit
                        break

                    # For every linker we add that is situated within the gap and not displacing an original residue,
                    # reduce the number we need to add to the following residues by 1
                    if lnum < interres - 1:
                        increase_orig_num = max(0, increase_orig_num - 1)
                    linker_id = base_id + lnum
                    linker_node = REBprocessor.Node(amino_acid=PDBtool.one_to_three(element),
                                                    residue_id=linker_id,
                                                    chain=new_nlist[-1].chain,
                                                    neighbor_prev=new_nlist[-1])
                    # Update previous node to point to linker
                    new_nlist[-1].neighbor_next = linker_node
                    # Add linker to list
                    new_nlist.append(linker_node)

                new_nlist[-1].neighbor_next = node
                node.neighbor_prev = new_nlist[-1]
                new_nlist.append(node)
            else:
                # Copy over original node with (potentially) adjusted residue id
                add_node = REBprocessor.Node(amino_acid=node.amino_acid,
                                             residue_id=node.residue_id + increase_orig_num,
                                             chain=node.chain,
                                             partners=node.partners,
                                             strength=node.strength,
                                             neighbor_prev=new_nlist[-1],
                                             orig_res_id=node.orig_res_id)
                new_nlist[-1].neighbor_next = add_node
                new_nlist.append(add_node)
        return new_nlist

    class SelectionStrategy(metaclass=ABCMeta):
        """
        Defines behavior a selection strategy needs to implement. Namely, given a list of nodes, select a subset as a
        peptide candidate.
        """

        def __init__(self):
            self.max_length = cm().get("peptide_generator.max_length")
            self.reb_energy_cutoff = cm().get("peptide_generator.reb_energy_cutoff")
            self.linker = cm().get("peptide_generator.linker")
            self.linking_force_length_limit = cm().get("peptide_generator.linking_force_length_limit")
            if self.reb_energy_cutoff > 0:
                logging.warning(f"You have set a positive energy cutoff ({self.reb_energy_cutoff}) "
                                f"- did you mean a negative value (attraction)?")
            self.length_damping = cm().get("peptide_generator.length_damping")
            self.ld_min_length = abs(cm().get("peptide_generator.length_damping_min_length"))
            self.ld_max_length = abs(cm().get("peptide_generator.length_damping_max_length"))
            self.ld_initial_mult = abs(cm().get("peptide_generator.length_damping_initial_mult"))
            self.ld_final_mult = abs(cm().get("peptide_generator.length_damping_final_mult"))
            self.ld_linear_stepping = cm().get("peptide_generator.length_damping_linear_stepping")
            if m := cm().get("peptide_generator.length_damping_mode"):
                if m.upper() == "QUADRATIC":
                    if self.ld_initial_mult > self.ld_final_mult:
                        self.damping_func = self._damping_factor_quadratic_penalty
                    else:
                        self.damping_func = self._damping_factor_quadratic_bonus
                elif m.upper() == "LINEAR":
                    self.damping_func = self._damping_factor_linear
                else:
                    logging.warning(f"Could not parse damping mode {m} for SelectionStrategy. "
                                    f"Using default (QUADRATIC)...")
                    self.damping_func = self._damping_factor_quadratic_penalty
            if self.ld_linear_stepping > 0:
                logging.warning(f"You have set a positive stepping ({self.ld_linear_stepping}) "
                                f"- did you perhaps mean a negative stepping (reduction of x% per step)?")
            if self.ld_initial_mult - (self.ld_max_length - self.ld_min_length) * self.ld_linear_stepping != \
                    self.ld_final_mult:
                logging.warning(f"Your configured settings for the linear damping factor yield a function that is not "
                                f"continuous. There will be a jump when reaching the length_damping_max_length!")

        @abstractmethod
        def reduce(self, from_chain: str, to_chain: str, nodes: List[REBprocessor.Node]) -> List[REBprocessor.Node]:
            raise NotImplementedError("Trying to use SelectionStrategy:reduce() from abstract base class!")

        def _damping_factor_quadratic_penalty(self, peptide_length: int) -> float:
            """
            Returns the quadratic damping factor.
            The damping factor is the value of a parabola open to the bottom, centered on the minimum length
            with a maximum of the initial bonus. Any length before the minimum length is equal to the initial
            multiplier. Any length after and including the max length is equal to the final multiplier.

                            ∧ (damping factor)
                            |
            initial bonus   |---*______
                            |          ‾‾‾‾----___
                            |                     ‾‾--_
            max penalty     |                          ‾*---------
                            |
            ----------------+---|-----------------------|----------> (peptide length)
                            |   min length               max length

            :param peptide_length: Length of the peptide, that the damping factor should be calculated for
            :return: Quadratic damping factor
            """
            return max(self.ld_final_mult,
                       (-1 * (max(peptide_length, self.ld_min_length) - self.ld_min_length) ** 2) *
                       max(0, (self.ld_initial_mult - self.ld_final_mult) /
                           (max(self.ld_min_length + 0.1, self.ld_max_length) - self.ld_min_length) ** 2) +
                       self.ld_initial_mult)

        def _damping_factor_quadratic_bonus(self, peptide_length: int) -> float:
            """
            Returns the quadratic damping factor.
            The damping factor is the value of a parabola centered on the minimum length with a maximum of the initial
            bonus. Any length before the minimum length is equal to the initial multiplier. Any length after and
            including the max length is equal to the final multiplier.

                            ∧ (damping factor)
                            |
            initial bonus   |                          _*---------
                            |                     __--‾
                            |          ____----‾‾‾
            max penalty     |---*‾‾‾‾‾‾
                            |
            ----------------+---|-----------------------|----------> (peptide length)
                            |   min length               max length

            :param peptide_length: Length of the peptide, that the damping factor should be calculated for
            :return: Quadratic damping factor
            """
            return max(self.ld_initial_mult,
                       (min(self.ld_final_mult, ((max(peptide_length, self.ld_min_length) - self.ld_min_length) ** 2) *
                            max(0.0, (self.ld_final_mult - self.ld_initial_mult) /
                                ((max(self.ld_min_length + 0.1, self.ld_max_length) - self.ld_min_length) ** 2)) +
                            self.ld_initial_mult)))

        def _damping_factor_linear(self, peptide_length: int) -> float:
            """
            Returns the linear damping factor.
            The damping factor is the initial bonus, discounted by the stepping for every residue above the min length,
            up to the max length.

                            ∧ (damping factor)
                            |
            initial bonus  -|---*___                (The slope equals the stepping)
                            |       ‾‾‾---___
                            |                ‾‾‾---___
            max penalty     |                         ‾‾‾*---------
                            |
            ----------------+---|------------------------|----------> (peptide length)
                            |   min length               max length


            :param peptide_length:
            :return: Linear damping factor
            """
            return (self.ld_initial_mult -
                    min((max(0, peptide_length - self.ld_min_length)), self.ld_max_length) * self.ld_linear_stepping)

    @final
    class GreedyExpand(SelectionStrategy):

        def __init__(self):
            super().__init__()
            self.always_include_direct_neighbors = cm().get(
                "peptide_generator.greedy_expand.always_include_direct_neighbors")
            self.max_side_extension = cm().get("peptide_generator.greedy_expand.max_side_extension")
            self.ignore_neighbors = cm().get("peptide_generator.greedy_expand.ignore_neighbors")
            self.custom_func = cm().get("peptide_generator.greedy_expand.custom_func")

        def reduce(self, from_chain: str, to_chain: str, nodes: List[REBprocessor.Node]) -> List[REBprocessor.Node]:
            if len(nodes) == 0:
                logging.error("Passed an empty list of nodes!")
                raise ValueError("Passed an empty list of nodes!")
            if any(_.chain not in from_chain for _ in nodes):
                logging.warning(f"There exist nodes in the passed nodes object that aren't on from_chain ({from_chain})!")

            def _include(add_to: list, n: REBprocessor.Node, curr_depth: int = 0) -> None:
                if self.custom_func:
                    import custom_funcs
                    return custom_funcs.greedy_expand_node_inclusion(self, add_to, n, curr_depth)
                if self.max_length == -1:
                    pass
                elif len(add_to) > self.max_length or n.strength_to(to_chain) >= 0:
                    return
                if self.always_include_direct_neighbors and curr_depth == 1:
                    add_to.append(n)
                damping_factor = self.damping_func(len(add_to)) if self.length_damping else 1
                if (n.strength_to(to_chain) * damping_factor < self.reb_energy_cutoff
                        and curr_depth <= self.max_side_extension
                        or (self.always_include_direct_neighbors and curr_depth == 1)):
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
                        if n.neighbor_prev.strength_to(to_chain) < n.neighbor_next.strength_to(to_chain):
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
            elligible_nodes = [n for n in nodes if n.chain in from_chain and n.strength_to(to_chain) < 0]
            if len(elligible_nodes) == 0:
                logging.warning(f"Couldn't determine any residues that interacted "
                                f"significantly according to the specifications. Exiting...")
                sys.exit(1)
            for residue in sorted(elligible_nodes, key=lambda node: node.strength_to(to_chain), reverse=False):
                _include(final_residues, residue)
            return _SelectionStrategies.add_linkers(
                sorted(final_residues, key=lambda node: node.residue_id, reverse=False),
                self.linking_force_length_limit)

    @final
    class FragmentJoiner(SelectionStrategy):

        def __init__(self):
            super().__init__()
            self.mode = ""
            if cm().get("peptide_generator.fragment_joiner.use_abs_increase"):
                self.mode += "A"  # Use absolute energy increase as criterion
            if cm().get("peptide_generator.fragment_joiner.use_rel_increase"):
                self.mode += "R"  # Use relative energy increase as criterion
            if len(self.mode) == 2 and cm().get("peptide_generator.fragment_joiner.mixed_mode_strict"):
                self.mode += "S"  # When using both, use strict mode (both need to increase)
            self.min_abs_increase = cm().get("peptide_generator.fragment_joiner.min_abs_increase")
            self.min_rel_increase = cm().get("peptide_generator.fragment_joiner.min_rel_increase")
            if self.min_rel_increase >= 1:
                logging.warning(f"You have set the minimum relative increase to something bigger than one - "
                                f"did you perhaps mean {self.min_rel_increase - 1} for "
                                f"a {(self.min_rel_increase - 1) * 100}% increase?")
            self.lookahead = cm().get("peptide_generator.fragment_joiner.lookahead")
            self.length_flexibility = cm().get("peptide_generator.fragment_joiner.length_flexibility")
            self.join_distance_penalty = cm().get("peptide_generator.fragment_joiner.join_distance_penalty")
            self.join_penalty_factor = cm().get("peptide_generator.fragment_joiner.join_penalty_factor", 0.05)
            self.fully_join_fragments = cm().get("peptide_generator.fragment_joiner.fully_join_fragments")
            self.penalize_lone_residues = cm().get("peptide_generator.fragment_joiner.penalize_lone_residues")
            self.lone_residue_penalty = cm().get("peptide_generator.fragment_joiner.lone_residue_penalty")
            self.pad_lone_residues = cm().get("peptide_generator.fragment_joiner.pad_lone_residues")
            self.lone_residue_pad_range = cm().get("peptide_generator.fragment_joiner.lone_residue_pad_range")
            self.old_frag_combiner = cm().get("peptide_generator.fragment_joiner.old_frag_combiner")
            self.linker_stretch_factor = cm().get("peptide_generator.fragment_joiner.linker_stretch_factor", 4.0)
            self.custom_func = cm().get("peptide_generator.fragment_joiner.custom_func")
            self.ref_relax = cm().get("ref_relax")

        @final
        @dataclass
        class _Fragment:
            n_ter_resid: int
            n_ter_cacoords: Tuple[float, float, float]
            c_ter_resid: int
            c_ter_cacoords: Tuple[float, float, float]
            nlist: List[REBprocessor.Node]
            score: float

            def __repr__(self):
                return f"[{' '.join([str(str(n) for n in self.nlist)])}]"

            def __str__(self):
                return self.__repr__()

            def __len__(self):
                return len(self.nlist)

        def reduce(self, from_chain: str, to_chain: str, nodes: List[REBprocessor.Node]) -> List[REBprocessor.Node]:
            if len(nodes) == 0:
                raise ValueError("Passed an empty list of nodes!")
            if any(_.chain not in from_chain for _ in nodes):
                logging.warning(
                    f"There exist nodes in the passed nodes object that aren't on from_chain ({from_chain})!")
            if self.ref_relax is None:
                if _ := cm().get("ref_relax"):
                    self.ref_relax = _
                else:
                    raise ValueError("Can't use FragmentJoiner strategy without the reference relaxed PDB!")

            def _include(n: REBprocessor.Node, to_chain: str, curr_strength: float, add_to_strength: float = 0,
                         curr_length: int = 0, ignore_cutoff: bool = False) -> bool:

                damping_factor = self.damping_func(curr_length) if self.length_damping else 1
                if self.custom_func:
                    import custom_funcs
                    return custom_funcs.fragment_joiner_node_inclusion(self, n, to_chain, curr_strength,
                                                                       add_to_strength,
                                                                       curr_length, ignore_cutoff, damping_factor)
                if not any(_ in to_chain for _ in nodes[residue_index].strength):  # Residue does not interact with target chain
                    return False
                abs_criterion = False
                rel_criterion = False
                node_strength = n.strength_to(to_chain)
                if node_strength == 0:  # Safety fallback, normally this shouldn't be encountered
                    node_strength = 0.00000001
                if curr_strength == 0:  # Safety fallback, normally this shouldn't be encountered
                    curr_strength = -0.00000001
                if "A" in self.mode and (ignore_cutoff or node_strength < self.reb_energy_cutoff):  # absolute increase
                    if add_to_strength + (node_strength * damping_factor) < self.min_abs_increase:
                        abs_criterion = True
                if "R" in self.mode and (ignore_cutoff or node_strength < self.reb_energy_cutoff):  # relative increase
                    new_strength = curr_strength + add_to_strength + (node_strength * damping_factor)
                    diff = new_strength - curr_strength
                    rel_change = diff / curr_strength
                    # We need to invert the % change if we're starting from a positive strength value,
                    # so that decreasing the strength value is good (positive percentage)
                    rel_change *= -1 if curr_strength > 0 else 1
                    if rel_change >= self.min_rel_increase:
                        rel_criterion = True
                return (("S" in self.mode and abs_criterion and rel_criterion) or  # Strict mode, both have to increase
                        ("S" not in self.mode and (abs_criterion or rel_criterion)))

            fragments = {}
            curr_fragment = []
            lookahead_buffer = []
            fragment_strength = -0.00000001  # Need to start with an epsilon, otherwise div/0 error in include

            # Forward scan
            for residue_index in range(len(nodes)):
                if nodes[residue_index] in curr_fragment:  # Already added to current fragment
                    continue

                added_residue = False
                lookahead_strength_buf = 0.0
                for i in range(self.lookahead + 1):
                    if nodes[residue_index + i].chain != from_chain:  # Only use residues that are not on target chain
                        break
                    if residue_index + i >= len(nodes):  # Don't lookahead past end of list
                        break
                    if len(curr_fragment) == 0 and i == 0 and nodes[residue_index + i].strength_to(to_chain) >= 0:
                        # Don't do lookahead if we are already starting from a node that doesn't interact with to_chain
                        break
                    # Do we add the current residue to our fragment?
                    if _include(nodes[residue_index + i], to_chain=to_chain, curr_strength=fragment_strength,
                                add_to_strength=lookahead_strength_buf,
                                curr_length=len(curr_fragment) + len(lookahead_buffer)):
                        added_residue = True
                        for r in lookahead_buffer:
                            curr_fragment.append(r)
                            fragment_strength += r.strength_to(to_chain, default=0)
                        curr_fragment.append(nodes[residue_index + i])
                        fragment_strength += nodes[residue_index + i].strength_to(to_chain, default=0)
                        break
                    else:  # Temporarily store current residue and look ahead
                        lookahead_strength_buf += nodes[residue_index + i].strength_to(to_chain, default=0)
                        lookahead_buffer.append(nodes[residue_index + i])
                lookahead_buffer = []
                if added_residue:  # During lookahead we added residue, we may proceed with current fragment
                    continue
                else:
                    if len(curr_fragment) > 0:
                        # We have a fragment and didn't find a suitable addition to the fragment during lookahead
                        if len(curr_fragment) == 1:  # fragment consists of only one residue
                            if self.pad_lone_residues:
                                neighbors = REBprocessor.Node.get_neighbors(curr_fragment[0],
                                                                            self.lone_residue_pad_range)
                                curr_fragment = []
                                fragment_strength = 0
                                for residue in neighbors:
                                    curr_fragment.append(residue)
                                    fragment_strength += residue.strength_to(to_chain, default=0)
                            if self.penalize_lone_residues:
                                fragment_strength += self.lone_residue_penalty

                        fragments[str(curr_fragment[0])] = \
                            (curr_fragment, fragment_strength)  # Save current fragment with aggregate energy
                        curr_fragment = []
                        fragment_strength = -0.00000001

            # FIXME: Alignment seems spurious, doesn't actually align any subseq of pep with its original position??

            # fragments have format:
            # {
            #   '1A': (
            #       [1A, 2A, 3A, 4A, ...],
            #       -0.85637828377...
            #   ),
            #   '19A': (
            #       [19A, 20A],
            #       -2.55555555...
            #   ),
            #   ...
            # }

            final_peptide = []

            if not self.old_frag_combiner:
                # Get coords of fragment N and C termini
                fraglist = []
                for _, fragment in fragments.items():
                    # PDB is N -> C
                    fraglist.append(_SelectionStrategies.FragmentJoiner._Fragment(
                        n_ter_resid=fragment[0][0].residue_id,
                        n_ter_cacoords=PDBtool.get_alphacarbon(self.ref_relax, fragment[0][0]),
                        c_ter_resid=fragment[0][-1].residue_id,
                        c_ter_cacoords=PDBtool.get_alphacarbon(self.ref_relax, fragment[0][-1]),
                        nlist=fragment[0],
                        score=fragment[1]
                    ))

                # get potential combinations
                def euclidean(coords1: Tuple[float, float, float], coords2: Tuple[float, float, float]) -> float:
                    return sqrt((coords1[0] - coords2[0]) ** 2 +
                                (coords1[1] - coords2[1]) ** 2 +
                                (coords1[2] - coords2[2]) ** 2)

                allowdist = len(self.linker) * self.linker_stretch_factor

                def grow(fragment_list: List, curr_combination: List, combination_list: List):
                    # TODO: Maybe use set instead of list for combination list so that duplicates get discarded?
                    # Get best possible n terminus extension of current fragment combination
                    n_partners = [f for f in fragment_list if f not in curr_combination and f != curr_combination[0]
                                  and euclidean(curr_combination[0].n_ter_cacoords, f.c_ter_cacoords) <= allowdist]
                    best_score = 10000000
                    best_n_frag = None
                    for _ in n_partners:
                        if _.score < best_score:
                            best_score = _.score
                            best_n_frag = _

                    # Get best possible c terminus extension of current fragment combination
                    c_partners = [f for f in fragment_list if f not in curr_combination and f != curr_combination[-1]
                                  and euclidean(curr_combination[-1].c_ter_cacoords, f.n_ter_cacoords) <= allowdist]
                    best_score = 10000000
                    best_c_frag = None
                    for _ in c_partners:
                        if _.score < best_score:
                            best_score = _.score
                            best_c_frag = _

                    # Identify whether to extend on N- or C-terminus side
                    if best_n_frag is not None and best_c_frag is not None:
                        use_frag = best_n_frag if best_n_frag.score < best_c_frag.score else best_c_frag
                    elif best_n_frag is None and best_c_frag is not None:
                        use_frag = best_c_frag
                    elif best_n_frag is not None and best_c_frag is None:
                        use_frag = best_n_frag
                    else:
                        fragstring = ""
                        for _ in curr_combination:
                            for __ in _.nlist:
                                fragstring = fragstring + str(__) + " "
                        fragstring = fragstring[:-1]  # remove trailing whitespace
                        logging.debug(f"Didn't find any eligible fragments to extend current fragment combination "
                                      f"({fragstring}).")
                        combination_list.append(curr_combination)
                        return

                    # Do we extend into the N-terminus region?
                    if use_frag.n_ter_resid < curr_combination[0].n_ter_resid:
                        # Would adding the current fragment exceed the maximum length?
                        if sum([len(_) for _ in curr_combination]) + len(use_frag) >= self.max_length:
                            # Yes, only use subset of residues until length limit is hit
                            use_num = self.max_length - sum([len(_) for _ in curr_combination])
                            subset = use_frag.nlist[-1 * use_num:]
                            use_frag = self._Fragment(
                                n_ter_resid=subset[0].residue_id,
                                n_ter_cacoords=PDBtool.get_alphacarbon(self.ref_relax, subset[0]),
                                c_ter_resid=subset[-1].residue_id,
                                c_ter_cacoords=PDBtool.get_alphacarbon(self.ref_relax, subset[-1]),
                                nlist=subset,
                                score=sum([_.strength_to(to_chain, default=0) for _ in subset])
                            )
                            curr_combination = [use_frag] + curr_combination
                            combination_list.append(curr_combination)
                            return
                        # There is still room to grow the current combination of fragments
                        curr_combination = [use_frag] + curr_combination
                        grow(fragment_list, curr_combination, combination_list)
                    else:
                        # Would adding the current fragment exceed the maximum length?
                        if sum([len(_) for _ in curr_combination]) + len(use_frag) >= self.max_length:
                            # Yes, only use subset of residues until length limit is hit
                            use_num = self.max_length - sum([len(_) for _ in curr_combination])
                            subset = use_frag.nlist[:use_num]
                            use_frag = self._Fragment(
                                n_ter_resid=subset[0].residue_id,
                                n_ter_cacoords=PDBtool.get_alphacarbon(self.ref_relax, subset[0]),
                                c_ter_resid=subset[-1].residue_id,
                                c_ter_cacoords=PDBtool.get_alphacarbon(self.ref_relax, subset[-1]),
                                nlist=subset,
                                score=sum([_.strength_to(to_chain, default=0) for _ in subset])
                            )
                            curr_combination.append(use_frag)
                            combination_list.append(curr_combination)
                            return
                        # There is still room to grow the current combination of fragments
                        curr_combination.append(use_frag)
                        grow(fragment_list, curr_combination, combination_list)

                # get feasible combinations of fragments
                buf = []
                for frag in fraglist:
                    grow(fragment_list=fraglist, curr_combination=[frag], combination_list=buf)
                if len(buf) == 0:
                    logging.error("Couldn't identify any fragment combinations! This is very likely an internal error, "
                                  "please create an issue on the VIPER GitHub.")
                    raise ValueError(
                        "Couldn't identify any fragment combinations! This is very likely an internal error, "
                        "please create an issue on the VIPER GitHub.")
                use_combination = buf[0]
                for comb in buf[1:]:
                    if sum([_.score for _ in comb]) < sum(_.score for _ in use_combination):
                        use_combination = comb
                final_peptide = []
                for frag in use_combination:
                    final_peptide += frag.nlist

            else:
                fragdict = sorted(fragments.items(), key=lambda _: _[1][1], reverse=False)
                # Apply inter-fragment distance penalty, if set
                if dist := self.join_distance_penalty:
                    pairterms = {}
                    # Get interfragment distances
                    for k1, fragment1 in fragdict:
                        pairterms[k1] = {}
                        for k2, fragment2 in fragdict:
                            if k2 == k1:
                                continue
                            pairterms[k1][k2] = \
                                PDBtool.get_dist_closest_atom(cm().get("ref_relax"), fragment1[0], fragment2[0])[0]
                    use_adjusted = []
                    # Actually applies distance penalties
                    # Assumes we're always including fragment no. 1 - TODO: don't assume this, test combinations
                    for k, frag in fragdict:
                        if k == fragdict[0][0]:
                            use_adjusted.append((frag[0][0], (frag[0], frag[1])))
                            continue
                        use_adjusted.append(
                            (frag[0][0],
                             (frag[0],
                              frag[1] * max(1 - max(pairterms[fragdict[0][0]][k] - dist, 0) * self.join_penalty_factor,
                                            0.01))))
                    fragdict = sorted(use_adjusted, key=lambda _: _[1][1], reverse=False)
                for k, fragment in fragdict:
                    stop = False
                    if len(final_peptide) + len(fragment[0]) + self.length_flexibility >= self.max_length:
                        continue
                    for residue in fragment[0]:
                        if not self.fully_join_fragments and len(
                                final_peptide) >= self.max_length + self.length_flexibility:
                            stop = True
                            break
                        else:
                            final_peptide.append(residue)
                    if stop:
                        break

            logging.debug(f"Final joined peptide fragments are: {pprint.pformat(final_peptide, compact=True)}")
            return _SelectionStrategies.add_linkers(
                sorted(final_peptide, key=lambda node: node.residue_id, reverse=False), self.linking_force_length_limit)

    @staticmethod
    def get_strategy(strategy: str = None) -> Type[_SelectionStrategies.SelectionStrategy]:
        strat = cm().get("peptide_generator.use_strategy")

        # a = test.get("peptide_generator.use_strategy")

        if strategy is not None:
            strat = strategy.upper()
        # strat = strategy.upper() if strategy is not None else cm.get("peptide_generator.use_strategy")
        if cm().get("peptide_generator.custom_strategy"):
            import custom_funcs
            return custom_funcs.CustomSelectionStrategy
        if strat == "GREEDYEXPAND" or strat == "GREEDY_EXPAND":
            return _SelectionStrategies.GreedyExpand
        elif strat == "FRAGMENTJOINER" or strat == "FRAGMENT_JOINER":
            return _SelectionStrategies.FragmentJoiner
        else:  # Default is FragmentJoiner strategy
            return _SelectionStrategies.FragmentJoiner
