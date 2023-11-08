from __future__ import annotations

import csv
import logging
import os.path
import re
import statistics
import subprocess
import sys
from dataclasses import dataclass, field
from enum import IntEnum
from io import StringIO
from pathlib import Path
from pprint import pformat
from typing import Union, List, final, Tuple, Dict

import pandas as pd

import ConfigManager

cm = ConfigManager.ConfigManager.get_cm


class RosettaWrapper:
    apps: dict = {}
    base_path: str = None

    def __init__(self):
        """
        Initialized the RosettaWrapper by autodiscovering all available Rosetta Apps in the given rosetta_config.path
        variable. If this variable is not set, raise a LookupError.
        Also sets up the folders to output files from Rosetta runs into.
        """
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
        else:
            raise LookupError("In order to use the RosettaWrapper, you must set the rosetta_config.path setting!")

        # Make directories
        self.base_path = os.path.join(cm().get("results_path"),
                                      os.path.normpath(cm().get("rosetta_config.path_out"))) + os.sep
        os.makedirs(self.base_path, exist_ok=True)
        os.makedirs(self.base_path + "run_configs", exist_ok=True)
        os.makedirs(self.base_path + "docking", exist_ok=True)
        os.makedirs(self.base_path + "prepack", exist_ok=True)
        os.makedirs(self.base_path + "refine", exist_ok=True)
        os.makedirs(self.base_path + "relax", exist_ok=True)

    def make_dir(self, p: Union[str, Path, List[str]]) -> str:
        """
        Creates a subdir in the folder for Rosetta output files.

        :param p: This may be either a path to directory or a name for a subdirectory, or a list of names of
            subdirectories to be nested. I.e. ["a", "b"] will create "basepath/a/b/"
        :return: The path that was created
        """
        use_path = self.base_path
        if isinstance(p, list):
            if len(p) == 0:
                raise ValueError("Can't pass empty list!")
            for subdir in p:
                use_path = os.path.join(use_path, subdir)
        else:
            use_path = os.path.join(use_path, p)
        os.makedirs(use_path, exist_ok=True)
        return use_path

    @staticmethod
    def make_options_file(filename: str, options: dict) -> str:
        """
        Writes an options dict to disk in the previously set up 'run_configs' folder as a Rosetta flag file.

        :param filename: Which filename to give the flag. Do NOT pass a full path to a file!
        :param options: The options to write to disk
            (Compare: https://www.rosettacommons.org/docs/latest/full-options-list)
        :return: The path to the written out flag file
        """
        use_path = os.path.join(cm().get("results_path"), os.path.normpath(cm().get("rosetta_config.path_out")),
                                "run_configs", filename)
        logging.info(f"Making run config {use_path}...")
        with open(use_path, "w") as out:
            for option, value in options.items():
                if value is not None:
                    use_value = str(value)
                    if isinstance(value, bool):
                        use_value = use_value.lower()
                    out.write(str(option) + " " + use_value + "\n")
                else:
                    out.write(str(option) + "\n")
        return use_path

    @staticmethod
    def _dispatch(app: str, flag: Union[str, Path]) -> None:
        """
        Runs a Rosetta app with a given flag.

        :param app: Which Rosetta app to use
        :param flag: Path to the flag file
        :return: None
        """
        flag = os.path.normpath(Path(flag))
        if ".mpi" in app:
            logging.debug(f"Running Rosetta with MPI...")
            res = subprocess.run(["mpirun", app, "@" + flag], capture_output=True, text=True)
            logging.debug("Stdout from Rosetta run: " + res.stdout)
            if len(res.stderr) == 0:
                if cm().get("permissive"):
                    logging.warning("The Rosetta run encountered errors! See: " + res.stderr)
                    logging.warning("Trying to continue anyway...")
                else:
                    raise RuntimeError("Rosetta threw errors: " + res.stderr)
        else:
            logging.debug(f"Running Rosetta with default...")
            res = subprocess.run([app, "@" + flag], capture_output=True, text=True)
            logging.debug("Stdout from Rosetta run: " + res.stdout)
            if len(res.stderr) == 0:
                if cm().get("permissive"):
                    logging.warning("The Rosetta run encountered errors! See: " + res.stderr)
                    logging.warning("Trying to continue anyway...")
                else:
                    raise RuntimeError("Rosetta threw errors: " + res.stderr)

    def run(self, flag: dict, options: Union[str, Path, dict], app: str = None) -> None:
        """
        This method provides a high level interface for running Rosetta apps.

        :param flag: A base flag dictionary, for a selection of a few sensible defaults, refer to the 'Flags' class
            at the end of this file
        :param options: Any (additional) options you want to provide. If they exist in 'flag', they will be overwritten.
            You may also pass a str/Path to an existing flag file. In this case you MUST pass a Rosetta app in 'app'
        :param app: Which app to use (optional). In case 'app' is specified in either 'flag' or 'options' has a value
            set for the key "app", that app will be used instead
        :return: None
        """
        if isinstance(options, str) or isinstance(options, Path):
            options = Path(options)
            if flag:
                logging.error("You cannot use a flag dictionary as well as a options file. Please supply a options"
                              "dict instead.")
                raise ValueError("You cannot use a flag dictionary as well as a options file. Please supply a options"
                                 "dict instead.")
            if not app:
                logging.error(f"When passing a path to an options file you also need to supply an "
                              f"application to run! Aborting...")
                raise ValueError("When passing a path to an options file you also need to supply an application "
                                 "to run! Aborting...")
            elif app not in self.apps.keys():
                logging.error(f"The application you wanted to run ({app}) could not be found! Aborting...")
                raise LookupError("The application you wanted to run ({app}) could not be found! Aborting...")
            if not os.path.isfile(os.path.normpath(options)):
                logging.error(f"The path '{options}' to the options file is not a file!")
            logging.info(f"Trying to run application '{app}' with options '{options}'...")
            RosettaWrapper._dispatch(self.apps[app], options)
        else:
            if flag:
                o = self.preprocess_options(flag, options)
            elif options:
                o = self.preprocess_options(options)
            else:
                logging.error("Need to set either flag or options!")
                raise ValueError("Need to set either flag or options!")
            use_app = app
            if a := o.pop("app", False):
                use_app = a
            use_app = self.apps.get(use_app, None)
            if use_app is None:
                raise ValueError("Couldn't determine which app to use! Please supply the name of the Rosetta app in"
                                 "either the flag dictionary or the options dictionary.")
            name = f"flag_{use_app}"
            flag_path = RosettaWrapper.make_options_file(filename=name, options=o)
            RosettaWrapper._dispatch(use_app, flag_path)

    def preprocess_options(self, flag: dict, override: dict = None, pdb: str = None) -> dict:
        """
        This method preprocesses a flag/options file for a run of Rosetta.
        First, all entries from override are copied over into flag, replacing existing values if necessary.
        Then the input file is set to 'pdb', or the PDB given at the start of VIPER as a fallback, if 'pdb' isn't set.

        Afterward, this method tries to fill in sensible values in 'flag' for certain keys and does some _very_
        rudimentary validation of the flag dictionary.

        :param flag: A base flag dictionary, for a selection of a few sensible defaults, refer to the 'Flags' class
            at the end of this file
        :param override: Any options in 'flag' that you wish to override
        :param pdb: Which pdb file to use as the input (optional)
        :return: A finalized options dictionary
        """
        logging.info("Preprocessing options file...")
        logging.debug("With flag:" + pformat(flag))
        logging.debug("With override:" + pformat(override))
        pdb_name = cm().get("PDB")
        if pdb:
            pdb_name = pdb
        if override:
            for k, v in override.items():
                flag[k] = v
        for k, v in flag.items():
            if flag["app"] in self.apps.keys():
                flag["app"] = self.apps[flag["app"]]
            else:
                logging.error(f"Could not preprocess options, because application '{flag['app']}' is unknown!")
                sys.exit(1)
            if "-in:file" in k and v is None:
                flag[k] = os.path.normpath(pdb_name)
                continue
            if "-out:" in k and v is None:
                if "file" in k and "fullatom" not in k:
                    if "relax" in flag["app"]:
                        flag[k] = os.path.join(self.base_path, "relax", pdb_name)
                    elif "prepack" in flag["app"]:
                        flag[k] = os.path.join(self.base_path, "prepack", pdb_name)
                    elif ("docking" in flag["app"] and not
                    ("-docking_local_refine" in flag.values() or "-docking_local_refine_min" in flag.values())):
                        flag[k] = os.path.join(self.base_path, "docking", pdb_name)
                    elif ("docking" in flag["app"] and
                          ("-docking_local_refine" in flag.values() or "-docking_local_refine_min" in flag.values())):
                        flag[k] = os.path.join(self.base_path, "refine", pdb_name)
                elif "path" in k:
                    if "relax" in flag["app"]:
                        flag[k] = os.path.join(self.base_path, "relax", pdb_name)
                    elif "prepack" in flag["app"]:
                        flag[k] = os.path.join(self.base_path, "prepack", pdb_name)
                    elif ("docking" in flag["app"] and not
                    ("-docking_local_refine" in flag.values() or "-docking_local_refine_min" in flag.values())):
                        flag[k] = os.path.join(self.base_path, "docking", pdb_name)
                    elif ("docking" in flag["app"] and
                          ("-docking_local_refine" in flag.values() or "-docking_local_refine_min" in flag.values())):
                        flag[k] = os.path.join(self.base_path, "refine", pdb_name)
                continue
            if "-run:constant_seed" in k:
                flag[k] = None
                flag["-run:jran"] = cm().get("rosetta_config.random_seed")
                continue
            if "-run:jran" in k:
                flag["-run:constant_seed"] = None
                flag[k] = cm().get("rosetta_config.random_seed")
                continue
            if "-nstruct" in k and v is None:
                # This one is dependent on specific context and can't be autodiscovered
                raise ValueError("Specified nstruct with no value, can't autodiscover appropriate value!")
            if "-backrub:ntrials" in k and v is None:
                flag[k] = cm().get("rosetta_config.relax_peptide_bb_ntrials")
            if "-dock_pert" in k and v is None:
                flag[k] = (f"{cm().get('rosetta_config.docking_translation')} "
                           f"{cm().get('rosetta_config.docking_rotation')}")
            if "-mh:path:scores_BB_BB" in k and v is None:
                flag[k] = os.path.join(os.path.normpath(cm().get("rosetta_config.path")), "main", "database",
                                       "additional_protocol_data", "motif_dock", "xh_16_")
        logging.info("Final flag: " + pformat(flag))
        return flag


class ScoreFileParser:
    """
    Parses score files from Rosetta apps
    """

    @staticmethod
    def read_in(score_file: Union[str, Path]) -> dict:
        """
        Reads a score file from a Rosetta run and returns its content as a dictionary of "pdb file": {values}

        :param score_file: Which score file to read
        :return: A dictionary containing the values for each pdb
        """
        logging.info(f"Reading in score file {score_file}...")
        formatted = ""
        # Reformat to CSV-like
        with open(score_file, "r") as i:
            pattern = re.compile(" +")
            for line in i:
                if "SEQUENCE:" in line:  # Skip first line
                    continue
                formatted += re.sub(pattern, ",", line)
        reader = csv.DictReader(StringIO(formatted))
        scores = {}
        try:
            for row in reader:
                scores[row["description"]] = {}
                for k, v in row.items():
                    if k in ["description", "SCORE:", ""]:  # Exclude fields with irrelevant data
                        continue
                    scores[row["description"]][k] = float(v)
        except csv.Error as e:
            if cm().get("permissive"):
                logging.warning(f"Encountered an error while processing score file {score_file}! Stacktrace: {e}")
            else:
                raise RuntimeError(e)
        return scores

    class ExtremumType(IntEnum):
        MINIMUM = 1
        MAXIMUM = -1

    @staticmethod
    def get_extremum(score_file: Union[str, Path], column: str,
                     extremum_type: ExtremumType = ExtremumType.MINIMUM) -> dict:
        """
        Gets the PDB that is associated with the specified extremum (minimum / maximum value) for the given column
        in the given score file.
        If no such PDB could be determined (usually if the specified colum can't be found in the score file), this will
        return None.

        :param score_file: Which Rosetta score file to analyze
        :param column: Which column to search for extremum
        :param extremum_type: Which extremum to search for (standard: ExtremumType.MINIMUM)
        :return: A dictionary in the format {"pdb_file": {dict of values for columns}} or None
        """
        scores = ScoreFileParser.read_in(score_file)
        curr_bound = 10000000000
        curr_best_pdb = None
        for pdb, values in scores.items():
            if v := values.get(column, None):
                if extremum_type * v < curr_bound:
                    curr_best_pdb = {pdb: scores[pdb]}
                    curr_bound = extremum_type * v
        return curr_best_pdb


class REBprocessor:
    """
    Parses output from Rosetta Residue Energy Breakdown
    """

    @staticmethod
    def read_in(breakdown_file: Union[str, Path]) -> Tuple[List, Dict]:
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
        """
        The Node dataclass holds all information from the residue energy breakdown parsing related to a specific
        residue.
        """
        amino_acid: str = None
        residue_id: int = None
        chain: str = None
        partners: List[tuple] = None  # (<partner residue id + chain id>, <interaction energy>)
        strength: dict = field(default_factory=lambda: ({}))
        neighbor_prev: REBprocessor.Node = None
        neighbor_next: REBprocessor.Node = None

        @staticmethod
        def get_neighbors(n: REBprocessor.Node, length: int, direction: int = 0, curr_depth: int = 0):
            """
            Returns the neighbors of the Node (residue) on the same chain, in order.
            For example, with this setup:

                A - B - C - D - E - F - G - H - J
                ---- increasing residue id ---->

            Calling this method on Node E with length 2 and direction 0 would give: [C, D, E, F, G].
            Calling this method on Node E with length 2 and direction < 0 would give: [C, D, E]

            :param n: The start node from which to gather the relative neighbors
            :param length: How many neighbors (inclusive) to get. I.e. "2" would give n + two neighbors
            :param direction: In which direction to look for neighbors. 0 looks in both directions (increasing and
                decreasing residue id), direction < 0 only considers neighbors with a smaller residue id,
                whereas direction > 0 only considers residues with a higher residue id
            :param curr_depth:
            :return:
            """
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

        def __contains__(self, item) -> bool:
            """
            Needed for custom 'in' behavior. Checks whether the residue id and chain id are the same as the input item.

            :param item: The item to compare against
            :return: True, if the concatenation of the own residue id and chain id is identical to item
            """
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


class Flags:
    residue_energy_breakdown: dict = {
        "app": "residue_energy_breakdown",
        "-in:file:s": None,
        "-out:file:silent": None,
        "-run:constant_seed": None,
        "-run:jran": None,
    }

    relax_complex_for_REB: dict = {
        "app": "relax",
        "in:file:s": None,
        "-out:path:all": None,
        "-out:suffix": "_relax_REB",
        "-run:constant_seed": None,
        "-run:jran": None,
        "-nstruct": 10,
        "-ex1": None,
        "-ex2": None,
        "-use_input_sc": None,
        "-no_his_his_pairE": None,
        "-no_optH": False,
        "-flip_HNQ": None,
        "-relax:thorough": None,
        "-relax:constrain_relax_to_start_coords": None,  # Pin backbone to known position
        "-relax:coord_constrain_sidechains": None,  # Pin heavy atoms in sidechains
        "-relax:ramp_constraints": False,  # Conserve constraints throughout run
    }

    relax_partner_protein: dict = {
        "app": "relax",
        "-in:file:s": None,
        "-out:path:all": None,
        "-out:suffix": "_relax_partner",
        "-run:constant_seed": None,
        "-run:jran": None,
        "-nstruct": None,  # relax_partner_runs
        "-ex1": None,
        "-ex2": None,
        "-use_input_sc": None,
        "-flip_HNQ": None,
        "-no_optH": None,
    }

    relax_peptide_normal_mode: dict = {
        "app": "rosetta_scripts",
        "-in:file:s": None,
        "-nstruct": None,  # relax_peptide_nmr_runs
        "-parser:protocol": os.path.join(os.path.realpath(__file__), "..", "..", "util", "rosetta_scripts",
                                         "normal_mode_relax.xml"),
        "-out:path:all": None,
        "-out:suffix": "_relax_peptide_normal_mode",
    }

    relax_peptide_backbone: dict = {
        "app": "relax",
        "-in:file:s": None,
        "-nstruct": None,  # relax_peptide_bb_runs
        "-backrub:ntrials": None,  # relax_peptide_bb_ntrials
        "-out:path:all": None,
        "-out:suffix": "_relax_peptide_backbone",
    }

    relax_peptide_fast: dict = {
        "app": "relax",
        "-in:file:s": None,
        "-nstruct": None,  # relax_peptide_fast_runs
        "-relax:thorough": None,
        "-out:path:all": None,
        "-out:suffix": "_relax_peptide_fast",
    }

    prepack_complex: dict = {
        "app": "docking_prepack_protocol",
        "-in:file:s": None,
        "-unboundrot": None,  # pdb
        "-nstruct": 1,
        "-partners": None,
        "-ensemble1": None,
        "-ensemble2": None,
        "-ex1": None,
        "-ex2aro": None,
        "-out:path:all": None,
        "-out:suffix": "_prepack_complex",
    }

    docking_ensemble: dict = {
        "app": "docking_protocol",
        "-in:file:s": None,
        "-unboundrot": None,  # pdb
        "-nstruct": None,  # docking_runs
        "-partners": None,
        "-ensemble1": None,
        "-ensemble2": None,
        "-dock_pert": None,  # docking_translation docking_rotation
        "-spin": None,
        "-detect_disulf": "true",
        "-rebuild_disulf": "true",
        "-ex1": None,
        "-ex2aro": None,
        "-docking_low_res_score": "motif_dock_score",
        "-mh:path:scores_BB_BB": None,
        "-mh:score:use_ss1": "false",
        "-mh:score:use_ss2": "false",
        "-mh:score:use_aa1": "true",
        "-mh:score:use_aa2": "true",
        "-out:path:all": None,
        "-out:suffix": "_docking_ensemble",
    }

    refine_local: dict = {
        "app": "docking_protocol",
        "-in:file:s": None,
        "-nstruct": None,  # refine_runs
        "-docking_local_refine": None,
        "-use_input_sc": None,
        "-ex1": None,
        "-ex2aro": None,
        "-out:file:fullatom": None,
        "-out:path:all": None,
        "-out:suffix": "_refine_local",
    }
