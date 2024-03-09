from __future__ import annotations

import collections
import copy
import csv
import logging
import os.path
import pprint
import re
import subprocess
import uuid
from dataclasses import dataclass, field
from enum import IntEnum
from io import StringIO
from pathlib import Path
from pprint import pformat
from typing import Union, List, final, Optional, Tuple

import pandas as pd

import ConfigManager
from util import Singleton

cm = ConfigManager.ConfigManager.get_instance


class RosettaWrapper(metaclass=Singleton._Singleton):
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
                elif ".default." in str(app.name) and len(self.apps.get(components[0], "")) == 0:
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
        logging.info(f"Making run config {use_path} ...")
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
            if "residue_energy_breakdown" in app:  # If REB is run with more than -np 1, there are duplicate entries
                res = subprocess.run(["mpirun", "-np", "1", app, "@" + flag], capture_output=True, text=True)
            else:
                res = subprocess.run(["mpirun", "-np", str(cm().get("rosetta_config.use_num_cores")), app, "@" + flag],
                                     capture_output=True, text=True)
            logging.debug("Stdout from Rosetta run: " + res.stdout)
            if len(res.stderr) > 0:
                if cm().get("permissive"):
                    logging.warning("The Rosetta run encountered errors! See: " + res.stderr)
                    logging.warning("Trying to continue anyway...")
                else:
                    raise RuntimeError("Rosetta threw errors: " + res.stderr)
        else:
            logging.debug(f"Running Rosetta with default...")
            res = subprocess.run([app, "@" + flag], capture_output=True, text=True)
            logging.debug("Stdout from Rosetta run: " + res.stdout)
            if len(res.stderr) > 0:
                if cm().get("permissive"):
                    logging.warning("The Rosetta run encountered errors! See: " + res.stderr)
                    logging.warning("Trying to continue anyway...")
                else:
                    raise RuntimeError("Rosetta threw errors: " + res.stderr)

    def run(self, flag: dict, options: Union[str, Path, dict], app: str = None, flag_suffix: str = None) -> None:
        """
        This method provides a high level interface for running Rosetta apps.

        :param flag: A base flag dictionary, for a selection of a few sensible defaults, refer to the 'Flags' class
            at the end of this file
        :param options: Any (additional) options you want to provide. If they exist in 'flag', they will be overwritten.
            You may also pass a str/Path to an existing flag file. In this case you MUST pass a Rosetta app in 'app'
        :param app: Which app to use (optional). In case 'app' is specified in either 'flag' or 'options' has a value
            set for the key "app", that app will be used instead
        :param flag_suffix: Which suffix to append to the resulting flag file. If None, a UUID will be used to prevent
            overwriting the configuration of a previous run
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
            app_name = None
            if flag:
                app_name = flag.get("app", None)
                if n := options.get("app", None):
                    app_name = n
                o = self.preprocess_options(flag, options)
            elif options:
                app_name = options.get("app", None)
                o = self.preprocess_options(options)
            else:
                logging.error("Need to set either flag or options!")
                raise ValueError("Need to set either flag or options!")
            use_app = app
            if a := o.pop("app", False):
                use_app = a
            if use_app is None:
                raise ValueError("Couldn't determine which app to use! Please supply the name of the Rosetta app in "
                                 "either the flag dictionary or the options dictionary.")
            if flag_suffix:
                name = f"flag_{app_name}_{flag_suffix}"
            else:
                name = f"flag_{app_name}_{uuid.uuid4()}"
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

        remove_parts = []
        if override:
            for k, v in override.items():
                flag[k] = v
        if flag["app"] in self.apps.keys():
            flag["app"] = self.apps[flag["app"]]
        else:
            logging.error(f"Could not preprocess options, because application '{flag['app']}' is unknown!")
            raise LookupError(f"Could not preprocess options, because application '{flag['app']}' is unknown!")
        for k, v in flag.items():
            # TODO: Maybe just check for absence of "in" and "out" params instead of trying to fill it in?
            if "-in:file" in k and v is None:
                if pdb:
                    flag[k] = os.path.normpath(pdb)
                else:
                    logging.error("No input file has been set!")
                    raise ValueError("No input file has been set!")
                continue
            elif "-in:file" in k and v is False:
                remove_parts.append(k)
            if "-out:" in k and v is None:
                if pdb is None:
                    logging.error("Can't derive path to output, because no input file has been set!")
                    raise ValueError("Can't derive path to output, because no input file has been set!")
                if "file" in k and "fullatom" not in k:
                    if "relax" in flag["app"]:
                        flag[k] = os.path.join(self.base_path, "relax", os.path.normpath(pdb))
                    elif "prepack" in flag["app"]:
                        flag[k] = os.path.join(self.base_path, "prepack", os.path.normpath(pdb))
                    elif ("docking" in flag["app"] and not
                    ("-docking_local_refine" in flag.values() or "-docking_local_refine_min" in flag.values())):
                        flag[k] = os.path.join(self.base_path, "docking", os.path.normpath(pdb))
                    elif ("docking" in flag["app"] and
                          ("-docking_local_refine" in flag.values() or "-docking_local_refine_min" in flag.values())):
                        flag[k] = os.path.join(self.base_path, "refine", os.path.normpath(pdb))
                elif "path" in k:
                    if "relax" in flag["app"]:
                        flag[k] = os.path.join(self.base_path, "relax", os.path.normpath(pdb))
                    elif "prepack" in flag["app"]:
                        flag[k] = os.path.join(self.base_path, "prepack", os.path.normpath(pdb))
                    elif ("docking" in flag["app"] and not
                    ("-docking_local_refine" in flag.values() or "-docking_local_refine_min" in flag.values())):
                        flag[k] = os.path.join(self.base_path, "docking", os.path.normpath(pdb))
                    elif ("docking" in flag["app"] and
                          ("-docking_local_refine" in flag.values() or "-docking_local_refine_min" in flag.values())):
                        flag[k] = os.path.join(self.base_path, "refine", os.path.normpath(pdb))
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
        for part in remove_parts:
            flag.pop(part)
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
        logging.debug(f"Reading in score file {score_file}...")
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
                    if k in ["description", "SCORE:", ""]:  # Only carry over actual score terms into the score dict
                        continue
                    scores[row["description"]][k] = float(v)
        except csv.Error as e:
            if cm().get("permissive"):
                logging.warning(f"Encountered an error while processing score file {score_file}! Stacktrace: {e}")
            else:
                raise RuntimeError(e)
        logging.debug(f"Scores are: {pprint.pformat(scores)}")
        return scores

    class ExtremumType(IntEnum):
        MINIMUM = 1
        MAXIMUM = -1

    @staticmethod
    def get_extremum(score_file: Union[str, Path], column: str,
                     extremum_type: ExtremumType = ExtremumType.MINIMUM) -> Optional[tuple]:
        """
        Gets the PDB that is associated with the specified extremum (minimum / maximum value) for the given column
        in the given score file.
        If no such PDB could be determined (usually if the specified colum can't be found in the score file), this will
        return None.

        :param score_file: Which Rosetta score file to analyze
        :param column: Which column to search for extremum
        :param extremum_type: Which extremum to search for (standard: ExtremumType.MINIMUM)
        :return: A tuple in the format ("path_to_pdb_file", {dict of values for columns}) or None
        """
        logging.debug(f"Getting extremum for file {score_file} ...")
        scores = ScoreFileParser.read_in(score_file)
        curr_bound = 10000000000
        curr_best_pdb = None
        for pdb, values in scores.items():
            if v := values.get(column, None):
                if extremum_type * v < curr_bound:
                    logging.debug(f"{pdb} is new best value with score {scores[pdb]}")
                    curr_best_pdb = (pdb, scores[pdb])
                    curr_bound = extremum_type * v
        logging.debug(f"Found extremum: {curr_best_pdb}")
        return curr_best_pdb


class REBprocessor:
    """
    Parses output from Rosetta Residue Energy Breakdown
    """

    @staticmethod
    def process_multipose(breakdown_file: Union[str, Path]) -> List[Node]:
        """
        Helper method to combine parsing of (multipose) residue energy breakdown file and calculate
        average energies per residue across all poses. Assumes that all poses are of the same molecules!

        :param breakdown_file: Path to breakdown file to be read
        :return: A list of Node objects with averaged interaction energies
        """
        try:
            breakdown_file = Path(breakdown_file).resolve(strict=True)
        except (FileNotFoundError, RuntimeError) as e:
            logging.error(f"Couldn't resolve a path you provided. Stacktrace: {e}")
            raise e
        return REBprocessor.get_avg_nodes(REBprocessor.read_in(breakdown_file))

    @staticmethod
    def read_in(breakdown_file: Union[str, Path]) -> collections.OrderedDict[List[Node]]:
        """
        Parses the output from a residue energy breakdown run and prepares it for further analysis.
        Can read in breakdown from file with multiple poses, returns them as a dict in the form of:
        {pose_id: [nodes...], ...}
        This assumes that every pose in the file is of the same molecule(s)!

        :param breakdown_file: A Rosetta Energy Breakdown score file
        :return: A dictionary of total Rosetta energies per residue per pose
        """
        logging.info(f"Reading in residue energy breakdown file {breakdown_file}...")
        formatted = ""
        # Reformat to CSV-like
        with open(breakdown_file, "r") as i:
            pattern = re.compile(" +")
            for line in i:
                formatted += re.sub(pattern, ",", line)
        breakdown_csv = pd.read_csv(StringIO(formatted),
                                    usecols=lambda c: c.upper() in ["POSE_ID", "PDBID1", "RESTYPE1", "PDBID2",
                                                                    "RESTYPE2",
                                                                    "TOTAL"])

        # Get amino acids on chains and interactions per amino acid
        pose_list = collections.OrderedDict()
        seen = {}
        for row in breakdown_csv.itertuples(index=False, name="Interaction"):
            if row.pose_id not in pose_list:
                pose_list[row.pose_id] = []
                seen[row.pose_id] = {}
            if row.restype2 == "onebody":  # Ignore onebody interactions
                continue
            # Create Node and parse data
            if row.pdbid1 not in seen[row.pose_id]:
                node = REBprocessor.Node(amino_acid=row.restype1, residue_id=int(row.pdbid1[:-1]), chain=row.pdbid1[-1],
                                         partners=[(row.pdbid2, float(row.total))],
                                         orig_res_id=(int(row.pdbid1[:-1]), row.pdbid1[-1]))
                if row.pdbid2[-1] in node.strength:
                    node.strength[row.pdbid2[-1]] += row.total
                else:
                    node.strength[row.pdbid2[-1]] = row.total
                seen[row.pose_id][row.pdbid1] = node
                pose_list[row.pose_id].append(node)
            else:  # Update Node with parsed data
                node = seen[row.pose_id][row.pdbid1]
                node.partners.append((row.pdbid2, float(row.total)))
                if row.pdbid2[-1] in node.strength:
                    node.strength[row.pdbid2[-1]] += row.total
                else:
                    node.strength[row.pdbid2[-1]] = row.total
            continue

        # Update neighbors
        for pose_id in seen:
            for residueid, node in seen[pose_id].items():
                prev_id = str(int(residueid[:-1]) - 1) + residueid[-1]
                if prev_id in seen[pose_id]:
                    node.neighbor_prev = seen[pose_id][prev_id]
                next_id = str(int(residueid[:-1]) + 1) + residueid[-1]
                if next_id in seen[pose_id]:
                    node.neighbor_next = seen[pose_id][next_id]

        return pose_list

    @staticmethod
    def get_avg_nodes(poses: collections.OrderedDict[List[Node]]) -> List[Node]:
        """
        Averages interaction energies for all residues in each pose and returns aggregated energies for each residue.

        :param poses: The output from the read_in() - A OrderedDict of {pose: List[Node], ...}
        :return: A list of nodes with averaged interaction energies
        """
        aggregate = []
        first = True
        num_poses = len(poses)
        if num_poses == 1:
            return poses.popitem()[1]
        for pose, nlist in poses.items():
            if first:
                # Get placeholder nodes for averaging energies
                aggregate = [REBprocessor.Node(n.amino_acid, n.residue_id, n.chain, n.partners.copy(),
                                               orig_res_id=copy.deepcopy(n.orig_res_id),
                                               strength=copy.deepcopy(n.strength)) for n in nlist]
                first = False
            else:
                # Sum all energies for each residue
                for n_index, node in enumerate(nlist):
                    for chain, energy in node.strength.items():
                        if chain in aggregate[n_index].strength:
                            aggregate[n_index].strength[chain] += energy
                        else:
                            aggregate[n_index].strength[chain] = energy
        # Divide by number of poses, save averaged energies from placeholder to actual nodes (which have neighbors etc.)
        use_pose = poses.popitem(last=False)[1]
        for n_index, node in enumerate(aggregate):
            for chain, energy in node.strength.items():
                use_pose[n_index].strength[chain] = node.strength[chain] / num_poses
        return use_pose

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
        orig_res_id: Tuple[int, str] = None

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

        def strength_to(self, to_chain: str, only_attr: bool = False, default: float = 10000000000) -> float:
            """
            Returns the aggregate strength for all interactions with all chains in to chain. If only_attr is True, only
            aggregate attractive interactions, if it is False (default), attractive and repulsive interactions get
            summed together and may cancel out.

            :param to_chain: A str of all concatened chain ids to aggregate the strength for, i. e. ACFH
            :param only_attr: Whether to only aggregate attractive interactions (default: false)
            :param default: What value to return if the interactions strength is 0, i. e. the residue doesn't interact
                with the chains in question
            :return: The aggregate interaction strength with the chains (or default)
            """
            if only_attr:
                interaction_strength = sum([self.strength.get(_, 0) for _ in to_chain if self.strength.get(_, 0) < 0])
            else:
                interaction_strength = sum([self.strength.get(_, 0) for _ in to_chain])
            return interaction_strength if interaction_strength != 0 else default

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


@final
class Flags(metaclass=Singleton._Singleton):
    @property
    def residue_energy_breakdown(self) -> dict:
        return {
            "app": "residue_energy_breakdown",
            "-score:weights": "ref2015",
            "-out:file:silent": None,
            "-run:constant_seed": None,
            "-run:jran": None,
            "-out:no_color": True,
        }

    @property
    def interface_analyzer(self) -> dict:
        return {
            "app": "InterfaceAnalyzer",
            "-score:weights": "ref2015",
            "-compute_packstat": True,
            "-packstat::oversample": 100,
            "-tracer_data_print": False,
            "-pack_input": True,
            "-pack_separated": True,
            "-add_regular_scores_to_scorefile": True,
            "-sasa_calculator_probe_radius": 1.4,
            "-pose_metrics::interface_cutoff": 8.0,
            "-use_input_sc": None,
            "-run:constant_seed": None,
            "-run:jran": None,
            "-out:no_color": True,
        }

    @property
    def relax_pinned_positions(self) -> dict:
        return {
            "app": "relax",
            "-score:weights": "ref2015",
            "-in:file:s": None,
            "-out:path:all": None,
            "-out:suffix": "_relax_pinned",
            "-run:constant_seed": None,
            "-run:jran": None,
            "-nstruct": 10,
            "-ex1": None,
            "-ex2": None,
            "-use_input_sc": None,
            "-no_his_his_pairE": None,
            "-no_optH": False,
            "-flip_HNQ": None,
            "-relax:fast": None,
            "-relax:constrain_relax_to_start_coords": None,  # Pin backbone to known position
            "-relax:coord_constrain_sidechains": None,  # Pin heavy atoms in sidechains
            "-relax:ramp_constraints": False,  # Conserve constraints throughout run
            "-out:no_color": True,
        }

    @property
    def relax_base(self) -> dict:
        return {
            "app": "relax",
            "-score:weights": "ref2015",
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
            "-out:no_color": True,
        }

    @property
    def relax_normal_mode(self) -> dict:
        return {
            "app": "rosetta_scripts",
            "-score:weights": "ref2015",
            "-in:file:s": None,
            "-nstruct": None,  # relax_peptide_nmr_runs
            "-parser:protocol": os.path.join(os.path.realpath(__file__), "..", "..", "util", "rosetta_scripts",
                                             "normal_mode_relax.xml"),
            "-out:path:all": None,
            "-out:suffix": "_relax_normal_mode",
            "-out:no_color": True,
        }

    @property
    def relax_backbone(self) -> dict:
        return {
            "app": "relax",
            "-score:weights": "ref2015",
            "-in:file:s": None,
            "-nstruct": None,  # relax_peptide_bb_runs
            "-backrub:ntrials": None,  # relax_peptide_bb_ntrials
            "-out:path:all": None,
            "-out:suffix": "_relax_backbone",
            "-out:no_color": True,
        }

    @property
    def relax_thorough(self) -> dict:
        return {
            "app": "relax",
            "-score:weights": "ref2015",
            "-in:file:s": None,
            "-nstruct": None,  # relax_peptide_fast_runs
            "-relax:thorough": None,
            "-out:path:all": None,
            "-out:suffix": "_relax_peptide_fast",
            "-out:no_color": True,
        }

    @property
    def prepack_complex(self) -> dict:
        return {
            "app": "docking_prepack_protocol",
            "-score:weights": "ref2015",
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
            "-out:no_color": True,
        }

    @property
    def docking_ensemble(self) -> dict:
        return {
            "app": "docking_protocol",
            "-score:weights": "ref2015",
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
            "-out:no_color": True,
        }

    @property
    def refine_local(self) -> dict:
        return {
            "app": "docking_protocol",
            "-score:weights": "ref2015",
            "-in:file:s": None,
            "-nstruct": None,  # refine_runs
            "-docking_local_refine": None,
            "-use_input_sc": None,
            "-ex1": None,
            "-ex2aro": None,
            "-out:file:fullatom": None,
            "-out:path:all": None,
            "-out:suffix": "_refine_local",
            "-out:no_color": True,
        }
