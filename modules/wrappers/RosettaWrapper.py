#  for constant randomness use '-run:constant_seed' and '-run:jran <seed>' (default: 1111111)
import logging
import os.path
import subprocess
import sys
from typing import Union, Tuple

from ConfigManager import ConfigManager
from VIPER import configmanager as cm


class RosettaWrapper:
    apps: dict = {}
    base_path: str = None

    def __init__(self, config: ConfigManager = None):

        # Initialize apps
        valid = False
        if r := cm.get("rosetta_config.path"):
            # TODO: Is this robust enough to determine that the path is pointing to the rosetta executables?
            if "/main/source/bin" not in r:
                app_path = os.path.join(os.path.normpath(r), "main", "source", "bin") + os.sep
            else:
                app_path = os.path.normpath(r) + os.sep
            for app in os.listdir(app_path):
                components = app.split(".")
                # For now, only accept mpi or default versions, and prefer mpi
                if ".mpi" in app:
                    valid = True  # We have found at least one executable, path is valid
                    self.apps[components[0]] = app
                    continue
                elif ".default" in app:
                    valid = True  # We have found at least one executable, path is valid
                    self.apps[components[0]] = app
        elif not valid:
            logging.critical(
                f"Path to Rosetta has been deemed invalid! Either it has not been set, or no apps could be found!")

        # Make directories
        self.base_path = os.path.join(cm.get("results_path"),
                                      os.path.normpath(cm.get("rosetta_config.path_out"))) + os.sep
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
        with open(os.path.join(cm.get("results_path"), os.path.normpath(cm.get("rosetta_config.path_out")),
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
                    f"When passing a path to run rosetta you also need to supply an application to run! Aborting... ")
                sys.exit(1)
            elif app not in self.apps.keys():
                logging.error(f"The application you wanted to run ({app}) could not be found! Aborting...")
                sys.exit(1)
            if os.path.isfile(os.path.normpath(options)):
                logging.info(f"Trying to run application '{app}' with options '{options}'...")
                RosettaWrapper._dispatch(self.apps[app], options)
        elif isinstance(options, dict):
            o = self.preprocess_options(options)

        pass

    def preprocess_options(self, flag: dict, override:dict = None, pdb: str = None) -> dict:
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
                    flag[k] = cm.get("PDB")
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
                flag[k] = cm.get("rosetta_config.random_seed")
                continue
            if "-nstruct" in k and v is not None:
                flag[k] = cm.get("rosetta_config.relax_protein_runs")
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
