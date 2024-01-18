from abc import ABCMeta
from typing import Tuple, Any


class OptimizationStrategy(metaclass=ABCMeta):

    def __init__(self):
        pass

    def run(self) -> Tuple[Any, float]:
        """
        This function will be run to determine the best result. The return value should be a tuple of the best
        individual and the associated score.

        :return: A tuple of (best individual, associated score)
        """
        raise NotImplementedError("Trying to use OptimizationStrategy:run() from abstract base class!")
