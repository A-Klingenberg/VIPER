from __future__ import annotations

import typing
from typing import List

if typing.TYPE_CHECKING:
    import modules.stages.PeptideGenerator
    import modules.stages.optimize.GAStrategy


def greedy_expand_node_inclusion(config: modules.stages.PeptideGenerator._SelectionStrategies.GreedyExpand,
                                 add_to: list, n: modules.stages.PeptideGenerator.REBprocessor.Node,
                                 depth: int) -> None:
    return NotImplemented


def fragment_joiner_node_inclusion(config: modules.stages.PeptideGenerator._SelectionStrategies.FragmentJoiner,
                                   n: modules.stages.PeptideGenerator.REBprocessor.Node, to_chain: str,
                                   curr_strength: float, add_to_strength: float = 0, curr_length: int = 0,
                                   ignore_cutoff: bool = True, damping_factor: float = 1) -> bool:
    return NotImplemented


class CustomSelectionStrategy(modules.stages.PeptideGenerator._SelectionStrategies.SelectionStrategy):

    def __init__(self):
        super().__init__()

    def reduce(self, from_chain: str, to_chain: str, nodes: List[modules.stages.PeptideGenerator.REBprocessor.Node]) -> \
            List[modules.stages.PeptideGenerator.REBprocessor.Node]:
        return NotImplemented


def addin_mutate(ga: modules.stages.optimize.GAStrategy.GAStrategy,
                 population: modules.stages.optimize.GAStrategy.Population) -> List:
    return NotImplemented


def custom_scii(scii: float) -> float:
    return NotImplemented
