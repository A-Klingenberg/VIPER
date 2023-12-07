from __future__ import annotations

import logging
import pprint
import sys
from abc import abstractmethod, ABCMeta
from typing import List, final, Type

import ConfigManager
from modules.wrappers.RosettaWrapper import REBprocessor
from util import PDBtool

cm = ConfigManager.ConfigManager.get_cm


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
        linker = cm().get("peptide_generator.linker")
        if linker is None:
            linker = "GG"  # default is polyglycine
        first = True
        new_list = []
        for node in nlist:
            if first:
                new_list.append(REBprocessor.Node(amino_acid=node.amino_acid,
                                                  residue_id=node.residue_id,
                                                  chain=node.chain,
                                                  orig_res_id=node.orig_res_id))
                first = False
                continue
            if abs(node.residue_id - new_list[-1].residue_id) > 1:  # Jump in the sequence, need linker
                base_id = new_list[-1].residue_id
                if respect_length_limit and len(new_list) + len(linker) >= cm().get("peptide_generator.max_length"):
                    raise IndexError("Cannot insert")
                for element in linker:  # Add linkers
                    base_id += 1
                    linker_node = REBprocessor.Node(amino_acid=PDBtool.one_to_three(element),
                                                    residue_id=base_id,
                                                    chain=new_list[-1].chain,
                                                    neighbor_prev=new_list[-1])
                    new_list[-1].neighbor_next = linker_node
                    new_list.append(linker_node)
                new_list[-1].neighbor_next = node
                node.neighbor_prev = new_list[-1]
                new_list.append(node)
            else:
                add_node = REBprocessor.Node(amino_acid=node.amino_acid,
                                             residue_id=node.residue_id,
                                             chain=node.chain,
                                             neighbor_prev=new_list[-1],
                                             orig_res_id=node.orig_res_id)
                new_list[-1].neighbor_next = add_node
                new_list.append(add_node)
        return new_list

    class SelectionStrategy(metaclass=ABCMeta):
        """
        Defines behavior a selection strategy needs to implement. Namely, given a list of nodes, generate a peptide
        candidate.
        """

        def __init__(self):
            self.max_length = cm().get("peptide_generator.max_length")
            self.reb_energy_cutoff = cm().get("peptide_generator.reb_energy_cutoff")
            self.linker = cm().get("peptide_generator.linker")
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
            if not from_chain == nodes[0].chain:
                raise ValueError(
                    "Chains have not been properly defined, list of nodes does not correspond to from_chain.")

            def _include(add_to: list, n: REBprocessor.Node, curr_depth: int = 0) -> None:
                if self.custom_func:
                    import custom_funcs
                    return custom_funcs.greedy_expand_node_inclusion(self, add_to, n, curr_depth)
                if self.max_length == -1:
                    pass
                elif len(add_to) > self.max_length or to_chain not in n.strength:
                    return
                if self.always_include_direct_neighbors and curr_depth == 1:
                    add_to.append(n)
                damping_factor = self.damping_func(len(add_to)) if self.length_damping else 1
                if (n.strength.get(to_chain, 10000000000) * damping_factor < self.reb_energy_cutoff
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
                        if (n.neighbor_prev.strength.get(to_chain, 10000000000) <
                                n.neighbor_next.strength.get(to_chain, 10000000000)):
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
            elligible_nodes = [n for n in nodes if n.chain == from_chain and to_chain in n.strength]
            if len(elligible_nodes) == 0:
                logging.warning(f"Couldn't determine any residues that interacted "
                                f"significantly according to the specifications. Exiting...")
                sys.exit(1)
            for residue in sorted(elligible_nodes, key=lambda node: node.strength[to_chain], reverse=False):
                _include(final_residues, residue)
            return _SelectionStrategies.add_linkers(
                sorted(final_residues, key=lambda node: node.residue_id, reverse=False))

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
            self.fully_join_fragments = cm().get("peptide_generator.fragment_joiner.fully_join_fragments")
            self.penalize_lone_residues = cm().get("peptide_generator.fragment_joiner.penalize_lone_residues")
            self.lone_residue_penalty = cm().get("peptide_generator.fragment_joiner.lone_residue_penalty")
            self.pad_lone_residues = cm().get("peptide_generator.fragment_joiner.pad_lone_residues")
            self.lone_residue_pad_range = cm().get("peptide_generator.fragment_joiner.lone_residue_pad_range")
            self.custom_func = cm().get("peptide_generator.fragment_joiner.custom_func")

        def reduce(self, from_chain: str, to_chain: str, nodes: List[REBprocessor.Node]) -> List[REBprocessor.Node]:
            if len(nodes) == 0:
                raise ValueError("Passed an empty list of nodes!")
            if not from_chain == nodes[0].chain:
                raise ValueError(
                    "Chains have not been properly defined, list of nodes does not correspond to from_chain.")

            def _include(n: REBprocessor.Node, to_chain: str, curr_strength: float, add_to_strength: float = 0,
                         curr_length: int = 0, ignore_cutoff: bool = False) -> bool:

                damping_factor = self.damping_func(curr_length) if self.length_damping else 1
                if self.custom_func:
                    import custom_funcs
                    return custom_funcs.fragment_joiner_node_inclusion(self, n, to_chain, curr_strength,
                                                                       add_to_strength,
                                                                       curr_length, ignore_cutoff, damping_factor)
                if to_chain not in nodes[residue_index].strength:  # Residue does not interact with target chain
                    return False
                abs_criterion = False
                rel_criterion = False
                node_strength = n.strength.get(to_chain, 10000000000)
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
            fragment_strength = -0.00000001  # Need to start with an epsilon, otherwise div/0_old error in include

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
                    # Do we add the current residue to our fragment?
                    if _include(nodes[residue_index + i], to_chain=to_chain, curr_strength=fragment_strength,
                                add_to_strength=lookahead_strength_buf,
                                curr_length=len(curr_fragment) + len(lookahead_buffer)):
                        added_residue = True
                        for r in lookahead_buffer:
                            curr_fragment.append(r)
                            fragment_strength += r.strength.get(to_chain, 0)
                        curr_fragment.append(nodes[residue_index + i])
                        fragment_strength += nodes[residue_index].strength.get(to_chain, 0)
                        break
                    else:  # Temporarily store current residue and look ahead
                        lookahead_strength_buf += nodes[residue_index + i].strength.get(to_chain, 0)
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
                                    fragment_strength += residue.strength.get(to_chain, 0)
                            if self.penalize_lone_residues:
                                fragment_strength += self.lone_residue_penalty

                        fragments[str(curr_fragment[0])] = \
                            (curr_fragment, fragment_strength)  # Save current fragment with aggregate energy
                        curr_fragment = []
                        fragment_strength = -0.00000001

            # TODO: Improve upon simple, greedy strategy?
            #  Knapsack problem, how to fit fragments with strength and length into max_length?
            #  If close to max, consider only adding subsequence of fragment? How to select?
            final_peptide = []

            def _get_strength(entry: tuple) -> float:
                return entry[1][1]

            temp = sorted(fragments.items(), key=_get_strength, reverse=False)

            for k, fragment in temp:
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

            logging.debug(f"Peptide fragments are: {pprint.pformat(temp)}")

            return _SelectionStrategies.add_linkers(
                sorted(final_peptide, key=lambda node: node.residue_id, reverse=False))

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
