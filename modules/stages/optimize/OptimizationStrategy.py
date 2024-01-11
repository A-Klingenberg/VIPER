from abc import ABCMeta, abstractmethod
from typing import List, Callable, Any


class OptimizationStrategy(metaclass=ABCMeta):
    mutators: List[Callable] = []

    def __init__(self):
        pass
