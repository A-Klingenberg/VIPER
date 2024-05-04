from abc import ABCMeta
from pathlib import Path
from typing import Union, List

from modules.wrappers import RosettaWrapper


class StructureProvider(metaclass=ABCMeta):

    def __init__(self):
        pass

    def get_structure(self, protein: Union[str, List[RosettaWrapper.REBprocessor.Node]]) -> Union[str, Path]:
        """
        This function will be run to get the tertiary structure of a protein / peptide. The result shall be either the
        PDB as a str, or a path to a saved PDB.

        :param protein: The protein to get the tertiary structure for. Can be either a str of single letter notation
            amino acids, or a list of REBprocessor.Node objects.
        :return: Either the content of the PDB as a str, or the path to the saved PDB file.
        """
        raise NotImplementedError("Trying to use StructureProvider:get_structure() from abstract base class!")
