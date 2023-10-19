from __future__ import annotations

#  for constant randomness use '-run:constant_seed' and '-run:jran <seed>' (default: 1111111)
import logging
import os.path
import re
import statistics
import subprocess
import sys
from dataclasses import dataclass, field
from io import StringIO
from pathlib import Path
from typing import Union, List, final

import pandas as pd

from ConfigManager import ConfigManager

cm = ConfigManager.get_cm


class RosettaWrapper:
    apps: dict = {}
    base_path: str = None

    def __init__(self, config: ConfigManager = None):

        # Initialize apps
        if r := cm().get("rosetta_config.path"):
            # TODO: This might be fairly slow. Defer loading to (non-blocking) background process and continue?
            logging.info(f"Loading Rosetta apps...")
            app_paths = [p.resolve() for p in Path(r).glob("**/*") if
                         (".mpi." in p.name or ".default." in p.name) and ("release" in p.name or "debug" in p.name)]
            for app in app_paths:
                components = app.name.split(".")
                # For now, only accept mpi or default versions, and prefer mpi
                if ".mpi." in str(app.name):
                    logging.debug(f"Found app '{components[0]}' in path '{str(app)}'!")
                    self.apps[components[0]] = str(app)
                elif ".default." in str(app.name):
                    logging.debug(f"Found app '{components[0]}' in path '{str(app)}'!")
                    self.apps[components[0]] = str(app)

        # Make directories
        self.base_path = os.path.join(cm().get("results_path"),
                                      os.path.normpath(cm().get("rosetta_config.path_out"))) + os.sep
        os.makedirs(self.base_path, exist_ok=True)
        os.makedirs(self.base_path + "run_configs", exist_ok=True)
        os.makedirs(self.base_path + "docking", exist_ok=True)
        os.makedirs(self.base_path + "prepack", exist_ok=True)
        os.makedirs(self.base_path + "refine", exist_ok=True)
        os.makedirs(self.base_path + "relax", exist_ok=True)

    @staticmethod
    def make_options_file(filename: str, options: dict) -> None:
        # Compare: https://www.rosettacommons.org/docs/latest/full-options-list
        logging.info(f"Making run config {filename}...")
        with open(os.path.join(cm().get("results_path"), os.path.normpath(cm().get("rosetta_config.path_out")),
                               "run_configs") + os.sep + filename, "w") as out:
            for option, value in options.items():
                out.write(str(option) + " " + str(value) + "\n")

    @staticmethod
    def _dispatch(app: str, flag: str) -> None:
        res = None
        if ".mpi" in app:
            logging.debug(f"Running Rosetta with MPI...")
            res = subprocess.run(["mpirun", app, "@" + flag], capture_output=True)
        else:
            logging.debug(f"Running Rosetta with default...")
            res = subprocess.run([app, "@" + flag], capture_output=True)

    def run(self, options: Union[str, dict], app: str = None) -> None:
        if isinstance(options, str):
            if not app:
                logging.error(
                    f"When passing a path to an options file you also need to supply an application to run! Aborting...")
                sys.exit(1)
            elif app not in self.apps.keys():
                logging.error(f"The application you wanted to run ({app}) could not be found! Aborting...")
                sys.exit(1)
            if not os.path.isfile(os.path.normpath(options)):
                logging.error(f"The path '{options}' to the options file is not a file!")
            logging.info(f"Trying to run application '{app}' with options '{options}'...")
            RosettaWrapper._dispatch(self.apps[app], options)
        elif isinstance(options, dict):
            o = self.preprocess_options(options)

        pass

    def preprocess_options(self, flag: dict, override: dict = None, pdb: str = None) -> dict:
        for k, v in flag.items():
            if flag["app"] in self.apps.keys():
                flag["app"] = self.apps[flag["app"]]
            else:
                logging.error(f"Could not preprocess options, because application '{flag['app']}' is unknown!")
                sys.exit(1)
            if "-in:file" in k:
                if pdb:
                    flag[k] = os.path.normpath(pdb)
                else:
                    flag[k] = cm().get("PDB")
                continue
            if "-out:file" in k and v is not None:
                if "relax" in flag["app"]:
                    flag[k] = os.path.join(self.base_path, "relax", flag[k])
                elif "prepack" in flag["app"] and v is not None:
                    flag[k] = os.path.join(self.base_path, "prepack", flag[k])
                elif ("docking" in flag["app"] and v is not None and
                      not ("-docking_local_refine" in flag.values() or "-docking_local_refine_min" in flag.values())):
                    flag[k] = os.path.join(self.base_path, "docking", flag[k])
                elif ("docking" in flag["app"] and v is not None
                      and ("-docking_local_refine" in flag.values() or "-docking_local_refine_min" in flag.values())):
                    flag[k] = os.path.join(self.base_path, "refine", flag[k])
                continue
            if "-run:constant_seed" in k and v is not None:
                flag[k] = cm().get("rosetta_config.random_seed")
                continue
            if "-nstruct" in k and v is not None:
                flag[k] = cm().get("rosetta_config.relax_protein_runs")
                continue
        return flag


class Flags:
    residue_energy_breakdown: dict = {
        "app": "residue_energy_breakdown",
        "-in:file:s": None,
        "-out:file:silent": None,
        "-run:constant_seed": None,
        "-run:jran": None,
    }

    relax_viral_surface_protein: dict = {
        "app": "relax",
        "-in:file:s": None,
        "-nstruct": None,

    }


class REBprocessor:
    """
    Parses results from Residue Energy Breakdown
    """

    @staticmethod
    def read_in(breakdown_file: str) -> tuple:
        """
        Parses the output from a residue energy breakdown run and prepares it for further analysis.

        :param breakdown_file:
        :return:
        """
        totals = []
        logging.info(f"Reading in residue energy breakdown file {breakdown_file}...")
        formatted = ""
        # Reformat to CSV-like
        with open(breakdown_file, "r") as i:
            pattern = re.compile(" +")
            for line in i:
                formatted += re.sub(pattern, ",", line)
        csv = pd.read_csv(StringIO(formatted),
                          usecols=lambda c: c.upper() in ["PDBID1", "RESTYPE1", "PDBID2", "RESTYPE2", "TOTAL"])

        # Get amino acids on chains and interactions per amino acid
        nlist = []
        seen = {}
        for row in csv.itertuples(index=False, name="Interaction"):
            if row.restype2 == "onebody":
                continue
            # Create Node and parse data
            if row.pdbid1 not in seen:
                node = REBprocessor.Node(amino_acid=row.restype1, residue_id=int(row.pdbid1[:-1]), chain=row.pdbid1[-1],
                                         partners=[(row.pdbid2, float(row.total))])
                if row.pdbid2[-1] in node.strength:
                    node.strength[row.pdbid2[-1]] += row.total
                else:
                    node.strength[row.pdbid2[-1]] = row.total
                seen[row.pdbid1] = node
                nlist.append(node)
                totals.append(row.total)
            else:  # Update Node with parsed data
                node = seen[row.pdbid1]
                node.partners.append((row.pdbid2, float(row.total)))
                totals.append(row.total)
                if row.pdbid2[-1] in node.strength:
                    node.strength[row.pdbid2[-1]] += row.total
                else:
                    node.strength[row.pdbid2[-1]] = row.total
            continue

        # Update neighbors
        for residueid, node in seen.items():
            prev_id = str(int(residueid[:-1]) - 1) + residueid[-1]
            if prev_id in seen:
                node.neighbor_prev = seen[prev_id]
            next_id = str(int(residueid[:-1]) + 1) + residueid[-1]
            if next_id in seen:
                node.neighbor_next = seen[next_id]

        stats = {"mean": statistics.mean(totals),
                 "stdev": statistics.stdev(totals),
                 "min": min(totals),
                 "max": max(totals)}
        return nlist, stats

    @final
    @dataclass
    class Node:
        amino_acid: str = None
        residue_id: int = None
        chain: str = None
        partners: List[tuple] = None  # (<partner residue id + chain id>, <interaction energy>)
        strength: dict = field(default_factory=lambda: ({}))
        neighbor_prev: REBprocessor.Node = None
        neighbor_next: REBprocessor.Node = None

        @staticmethod
        def get_neighbors(n: REBprocessor.Node, length: int, direction: int = 0, curr_depth: int = 0):
            if curr_depth > length or n is None:
                return None
            if direction < 0:
                if p := REBprocessor.Node.get_neighbors(n.neighbor_prev, length, direction, curr_depth + 1):
                    return p + [n]
                else:
                    return [n]
            if direction > 0:
                if nex := REBprocessor.Node.get_neighbors(n.neighbor_next, length, direction, curr_depth + 1):
                    return [n] + nex
                else:
                    return [n]
            if direction == 0:  # Initialize
                buf = [n]
                if p := REBprocessor.Node.get_neighbors(n.neighbor_prev, length, -1, curr_depth + 1):
                    buf = p + buf
                if n := REBprocessor.Node.get_neighbors(n.neighbor_next, length, 1, curr_depth + 1):
                    buf = buf + n
                return buf

        def __contains__(self, item):
            return item == str(self.residue_id) + self.chain

        def __repr__(self):
            return str(self.residue_id) + self.chain + " " + str(self.amino_acid)

        def __int__(self):
            """
            This returns the residue id!

            :return: residue id
            """
            return self.residue_id

        def __str__(self):
            return str(self.residue_id) + self.chain
