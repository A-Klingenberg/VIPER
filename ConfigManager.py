"""
The ConfigManager parses, stores, and makes all configuration data available within VIPER.
"""
import argparse
import json
import logging
import os.path
import random
import sys
from pathlib import Path
from pprint import pformat
from typing import Optional, Any, Union

import numpy.random

import version
from util import Singleton


class ConfigManager(metaclass=Singleton._Singleton):
    command_args: dict = None
    file_config: dict = None
    program_path: str = None
    results_path: str = None
    log_path: str = None
    summary_logger: logging.Logger = None
    ref_relax: Union[Path, str] = None

    def __init__(self, args: argparse.Namespace = None,
                 base_path: str = os.path.abspath(os.path.join(os.path.dirname(__file__), f"..{os.sep}"))):
        """
        This method initializes the config manager for VIPER. It saves any arguments given to the program and sets any
        further settings, if necessary.
        Important: command-line settings take precedence over config.json settings!

        :param args: The command-line arguments given to VIPER, as parsed by argparse.ArgumentParser().parse_args()
        :param base_path: The base path of VIPER. Based on this, subfolders will be created and files will be saved.
        """
        self.program_path = os.path.normpath(base_path) + os.sep
        c_path = None
        if args:
            self.command_args = vars(args)
            if self.command_args["config"]:
                c_path = self.command_args["config"]
        if not c_path:
            # Fallback: If path to config file is not set via command line, expect config.json in same folder as main.py
            c_path = self.program_path + "config.json"

        try:
            with open(c_path, 'r') as config_file:
                conf_in: dict = json.load(config_file)
                # TODO: Is there a need to validate file_config?

                # recursively flatten JSON structure
                json_conf = {}

                # TODO: Support list entries?
                def flatten_level(d: dict, add_to: dict, prefix: str = None) -> None:
                    for k, v in d.items():
                        if isinstance(v, dict):
                            if prefix:
                                flatten_level(v, add_to, prefix=prefix + f"{str(k)}.")
                            else:
                                flatten_level(v, add_to, prefix=f"{str(k)}.")
                        else:
                            if prefix:
                                add_to[f"{prefix}{str(k)}"] = v
                            else:
                                add_to[str(k)] = v

                flatten_level(conf_in, json_conf)

                self.file_config = json_conf
        except OSError:
            print(f"Could not open config file at '{c_path}'.")
            if self.command_args and "permissive" not in self.command_args:
                print(f"Aborting...")
                sys.exit(1)
            else:
                # Note: Permissive is true by default
                print(f"Trying to continue with just command line arguments - this might cause errors down the line!")

        # Initialize main logger
        if lp := self.get("log_path"):
            self.log_path = os.path.normpath(self.program_path + lp) + os.sep
        else:  # Fallback
            self.log_path = self.program_path + "Logs" + os.sep

        # If the folder where the logs shall be saved doesn't exist, create it
        os.makedirs(self.log_path, exist_ok=True)

        # Summary
        summary_formatter = logging.Formatter("%(asctime)s | (%(module)s): %(message)s")
        self.summary_logger = logging.getLogger("summary")
        self.summary_logger.setLevel(logging.INFO)
        fhnd = logging.FileHandler(os.path.abspath(self.log_path + "summary.txt"))
        fhnd.setFormatter(summary_formatter)
        self.summary_logger.addHandler(fhnd)

        # Regular
        if self.get("verbose"):
            logging.basicConfig(level=logging.DEBUG,
                                filename=os.path.abspath(self.log_path + "log.txt"),
                                format="%(asctime)s [%(levelname)s] (%(filename)s:%(module)s): %(message)s",
                                force=True)
        else:
            logging.basicConfig(level=logging.INFO,
                                filename=os.path.abspath(self.log_path + "log.txt"),
                                format="%(asctime)s [%(levelname)s] (%(filename)s:%(module)s): %(message)s",
                                force=True)
        self.log_config()

        # Set output path
        p = None
        if (self.file_config and "results_path" in self.file_config.keys() and
                self.file_config["results_path"] is not None):
            p = self.file_config["results_path"]
        if (self.command_args and "results_path" in self.command_args.keys() and
                self.command_args["results_path"] is not None):
            p = self.command_args["results_path"]
        if p is None:
            p = "output"  # Fallback
        self.results_path = os.path.join(self.program_path, p) + os.sep

        # Seed randomness
        random.seed(self.get("random_seed"))
        numpy.random.seed(self.get("random_seed"))

        # Validate Rosetta path
        if r := self.get("rosetta_config.path"):
            # Try to make sure that Rosetta executables are findable in given directory
            apps = [p.resolve() for p in Path(r).glob("**/*") if
                    (".mpi." in p.name or ".default." in p.name) and ("release" in p.name or "debug" in p.name)]
            if len(apps) == 0:
                if self.get("permissive"):
                    logging.warning(f"Could not find any Rosetta executables in the path you provided! ({r}) "
                                    f"If you continue anyway and try to use Rosetta apps, you will run into errors!")
                else:
                    logging.error(
                        f"Could not find any Rosetta executables in the path you provided! ({r}) - Aborting...")
                    sys.exit(1)

    def get(self, setting: str, default: Any = None) -> Optional[Any]:
        """
        Gets a given setting, if it is set, otherwise returns None. For a list of all settings,
        see accompanying documentation.
        Important: command-line settings take precedence over config.json settings!

        :param setting: A string representing the name of the setting to be retrieved
        :param default: What to return if the setting cannot be retrieved
        :return: The value of the setting, or default (None) if that setting hasn't been configured
        """
        if setting == "program_path":
            return self.program_path
        if setting == "results_path":
            return self.results_path
        if setting == "ref_relax":
            return self.ref_relax
        if self.command_args and setting in self.command_args and self.command_args[setting] is not None:
            if isinstance(self.command_args[setting], str) and len(self.command_args[setting]) == 0:
                if default:
                    return default
                else:
                    return None
            return self.command_args[setting]
        elif self.file_config and setting in self.file_config and self.file_config[setting] is not None:
            if isinstance(self.file_config[setting], str) and len(self.file_config[setting]) == 0:
                if default:
                    return default
                else:
                    return None
            return self.file_config[setting]
        else:
            return default

    def log_config(self) -> None:
        """
        Logs all configuration settings to the log file.

        :return: None
        """
        self.summary_logger.info(f"Running VIPER version {version.VERSION}.")
        logging.info(f"Running Python ver. {sys.version}")
        logging.info(f"Running VIPER version {version.VERSION}.")
        logging.info(f"Using the following configuration:")
        if self.file_config:
            logging.info(f"---------- [ CONFIG FILE ] ----------")
            logging.info(pformat(self.file_config))
        if self.command_args:
            logging.info(f"---------- [ COMMAND LINE ARGUMENTS ] ----------")
            logging.info(pformat(self.command_args))
