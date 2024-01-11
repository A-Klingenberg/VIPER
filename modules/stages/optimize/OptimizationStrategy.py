from abc import ABCMeta
from typing import List, Callable


class OptimizationStrategy(metaclass=ABCMeta):
    mutators: List[Callable] = []

    def __init__(self):
        pass
