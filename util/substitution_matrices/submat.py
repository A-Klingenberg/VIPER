import json
import logging
import os.path
import pprint
import sys
from pathlib import Path
from typing import Union

import ConfigManager
from util import PDBtool

cm = ConfigManager.ConfigManager.get_instance


class SubMat:
    mat: dict = None

    def __init__(self, name: Union[str, Path], fallback="BLOSUM62_shifted.json"):
        use_fallback = False
        use_name = None
        if isinstance(name, Path):
            try:
                use_name = Path(name).resolve(strict=True)
            except (FileNotFoundError, RuntimeError) as e:
                logging.error(f"Could not resolve the path provided ({name}): {e}")
                use_fallback = True
        else:
            if not name.endswith(".json"):
                name = name + ".json"
            try:
                use_name = Path(name).resolve(strict=True)
            except (FileNotFoundError, RuntimeError) as e:
                logging.warning(f"The passed name ({name}) couldn't be parsed as a Path: {e}")
                logging.warning("Trying to find files with that name in substitution_matrices/ instead...")
            p = Path(os.path.join(cm().get("program_path"), "util", "substitution_matrices")).resolve(strict=True)
            matches = [_ for _ in p.glob(f"**/{name}")]
            if len(matches) == 0:
                logging.error(f"Could not find the substitution matrix you indicated ({name})")
                use_fallback = True
            else:
                logging.debug(f"Identified the following matching files: {pprint.pformat(matches, compact=True)}")
                for m in matches:
                    if m.name == name:  # exact match
                        use_name = m
                use_name = matches[0] if use_name is None else use_name  # No exact match found, use first match
                logging.debug(f"Identified substitution matrix {use_name}")
        if not use_fallback:
            try:
                with open(use_name, "r") as submat:
                    self.mat = json.load(submat)
            except (FileNotFoundError, IsADirectoryError, json.JSONDecodeError) as e:
                logging.error(f"Couldn't load substitution matrix {use_name}: {e}")
                sys.exit(1)
        else:
            fallback = os.path.join(cm().get("program_path"), "util", "substitution_matrices") + os.sep + fallback
            try:
                fallback = Path(fallback).resolve(strict=True)
                logging.warning(f"Using fallback: {fallback}")
                with open(fallback, "r") as submat:
                    self.mat = json.load(submat)
            except (FileNotFoundError, IsADirectoryError, RuntimeError, json.JSONDecodeError) as ef:
                logging.error(f"Couldn't load fallback matrix: {ef}")
                sys.exit(1)

    def get(self, key: str):
        """
        Gets the substitution biases for the given amino acids. No guarantee is made regarding the format of the result.
        Usually this should be a dictionary of {amino acid: numeric value}.

        :param key: For which amino acid to get the substitution values
        :return: Usually a dictionary of format {amino acid: numeric value}
        """
        if _ := self.mat.get(key):
            return _
        elif _ := self.mat.get(PDBtool.one_to_three(key)):
            return _
        elif _ := self.mat.get(PDBtool.three_to_one(key)):
            return _


class BLOSUM45(SubMat):
    # Taken from:
    # https://web.archive.org/web/20240102153211/https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM45
    # Entries for the BLOSUM45 matrix at a scale of ln(2)/3.0.
    def __init__(self):
        super().__init__("BLOSUM45.json")


class BLOSUM45_s(SubMat):
    def __init__(self):
        super().__init__("BLOSUM45_shifted.json")


class BLOSUM50(SubMat):
    # Taken from:
    # https://web.archive.org/web/20220803023229/https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM50
    # Entries for the BLOSUM50 matrix at a scale of ln(2)/3.0.
    def __init__(self):
        super().__init__("BLOSUM50.json")


class BLOSUM50_s(SubMat):
    def __init__(self):
        super().__init__("BLOSUM50_shifted.json")


class BLOSUM62(SubMat):
    # Taken from:
    # https://web.archive.org/web/20231122014358/https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM62
    # Entries for the BLOSUM62 matrix at a scale of ln(2)/2.0.
    def __init__(self):
        super().__init__("BLOSUM62.json")


class BLOSUM62_s(SubMat):
    def __init__(self):
        super().__init__("BLOSUM62_shifted.json")


class BLOSUM80(SubMat):
    # Taken from:
    # https://web.archive.org/web/20231118164317/http://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM80
    # Entries for the BLOSUM80 matrix at a scale of ln(2)/2.0.
    def __init__(self):
        super().__init__("BLOSUM80.json")


class BLOSUM80_s(SubMat):
    def __init__(self):
        super().__init__("BLOSUM80_shifted.json")


class BLOSUM90(SubMat):
    # Taken from:
    # https://web.archive.org/web/20220707144125/http://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM90
    # Entries for the BLOSUM90 matrix at a scale of ln(2)/2.0.
    def __init__(self):
        super().__init__("BLOSUM90.json")


class BLOSUM90_s(SubMat):
    def __init__(self):
        super().__init__("BLOSUM90_shifted.json")


class UNIFORM(SubMat):
    def __init__(self):
        super().__init__("UNIFORM.json")
