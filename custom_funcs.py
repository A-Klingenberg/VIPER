from __future__ import annotations

import operator
import typing
from typing import List, Union

if typing.TYPE_CHECKING:
    import modules.stages.PeptideGenerator
    import modules.stages.optimize.GAStrategy


def greedy_expand_node_inclusion(config: modules.stages.PeptideGenerator._SelectionStrategies.GreedyExpand,
                                 add_to: list, n: modules.stages.PeptideGenerator.REBprocessor.Node,
                                 depth: int) -> None:
    """
    Determines whether a node (residue) should be included when growing a fragment during the GreedyExpand selection
    strategy. 'add_to' is the current fragment, 'n' the residue in question, and 'depth' the current distance to the
    starting point of this fragment.

    :param config: The GreedyExpand instance that VIPER created. You can access settings like 'max_side_extension' via
        this object.
    :param add_to: The current fragment, with all residues added so far.
    :param n: The residue that is being considered for inclusion.
    :param depth: The distance in number of residues to the residue that started this fragment.
    :return: Nothing. Modify add_to as needed (by recursively calling this function), and return None when done.
    """
    return NotImplemented


def fragment_joiner_node_inclusion(config: modules.stages.PeptideGenerator._SelectionStrategies.FragmentJoiner,
                                   n: modules.stages.PeptideGenerator.REBprocessor.Node, to_chain: str,
                                   curr_strength: float, add_to_strength: float = 0, curr_length: int = 0,
                                   ignore_cutoff: bool = True, damping_factor: float = 1) -> bool:
    """
    Determines whether a node (residue) should be included when growing a fragment during the FragmentJoiner selection
    strategy. 'n' is the residue in question, 'to_chain' describes the VSP chain(s), 'curr_strength' describes the
    strength of the current fragment, 'add_to_strength' is a temporary buffer that holds the aggregte score for the
    lookahead window, 'curr_length' is the current length of the fragment in question, 'ignore_cutoff' describes whether
    the residue energy should be ignored, and 'damping_factor' is the precomputed damping factor that VIPER passes
    along to you.

    :param config: The FragmentJoiner instance that VIPER created. You can access settings like 'lookahead' via this
        object.
    :param n: The residue that is being considered for inclusion.
    :param to_chain: The VSP chain(s)
    :param curr_strength: The strength of the current fragment.
    :param add_to_strength: The aggregate strength of the lookahead window so far.
    :param curr_length: The length of the current fragment.
    :param ignore_cutoff: Whether to ignore the residue energy during inclusion.
    :param damping_factor: The damping factor to apply to this residue's score. VIPER precomputes this for you.
    :return: True or False, depending on whether the residue should be included or not.
    """
    return NotImplemented


class CustomSelectionStrategy(modules.stages.PeptideGenerator._SelectionStrategies.SelectionStrategy):
    """
    You can implement your own selection strategy here. Have a look at the selection strategies in PeptideGenerator to
    understand how these should be architected. Make sure your selection strategy implements the 'reduce' function below,
    as this is the interface selection strategies must implement for VIPER to work correctly with them.
    """

    def __init__(self):
        super().__init__()

    def reduce(self, from_chain: str, to_chain: str, nodes: List[modules.stages.PeptideGenerator.REBprocessor.Node]) -> \
            List[modules.stages.PeptideGenerator.REBprocessor.Node]:
        """
        Reduce a list of nodes (residues) to a candidate peptide (also a list of nodes).

        :param from_chain: The partner chain(s), from where residues should be taken.
        :param to_chain: The VSP chain(s).
        :param nodes: The list of nodes from the residue energy breakdown of a PDB of the complex.
        :return: A list of nodes (residues) representing the candidate peptide.
        """
        return NotImplemented


def addin_mutate(ga: modules.stages.optimize.GAStrategy.GAStrategy,
                 population: modules.stages.optimize.GAStrategy.Population) -> List:
    """
    This function can be called after the mutation step. Here you can directly modify each population of the generation.
    This function will be called separately for each population in the generation.

    :param ga: The GAStrategy instance that VIPER created. You can access settings like 'selection_mode' via this
        object's config dict.
    :param population: The population in question.
    :return: A list of individuals that should be added to the current population.
    """
    return NotImplemented


def custom_scii_bonus(scii: float) -> Union[
    Union[operator.mul, operator.add, operator.sub, operator.floordiv, operator.pow, operator.truediv], float]:
    """
    This function defines how the sSCII value is combined with the Rosetta dG_separated term.

    :param scii: The scii value for this candidate.
    :return: A tuple of a supported mathematical operator (such mul for multiplication) and the other operand that
        should be used when combining the scores. For example, a tuple of (operator.mul, 0.95) would mean that the
        Rosetta score will be multiplied with 0.95, so a 5% penalty.
    """
    return NotImplemented
