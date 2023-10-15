from typing import List

import modules.stages.PeptideGenerator


def greedy_expand_node_inclusion(config: modules.stages.PeptideGenerator._SelectionStrategies.GreedyExpand,
                                 add_to: list, n: modules.stages.PeptideGenerator.REBprocessor.Node,
                                 depth: int) -> None:
    return NotImplemented


def fragment_joiner_node_inclusion(config: modules.stages.PeptideGenerator._SelectionStrategies.FragmentJoiner,
                                   n: modules.stages.PeptideGenerator.REBprocessor.Node, curr_strength: float,
                                   add_to_strength: float = 0, curr_length: int = 0,
                                   ignore_cutoff: bool = True, damping_factor: float = 1) -> bool:
    return NotImplemented


class CustomSelectionStrategy(modules.stages.PeptideGenerator._SelectionStrategies.SelectionStrategy):

    def __init__(self):
        super().__init__()

    def reduce(self, from_chain: str, to_chain: str, nodes: List[modules.stages.PeptideGenerator.REBprocessor.Node]) -> \
            List[modules.stages.PeptideGenerator.REBprocessor.Node]:
        return NotImplemented
