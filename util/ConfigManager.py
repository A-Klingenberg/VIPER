"""
The ConfigManager parses, stores, and makes all configuration data available within VIPER.

Version: 2023-09-27
"""
import argparse
import json
import logging
import os.path
import sys
from pprint import pformat
from typing import Optional, Any

import VIPER


class ConfigManager:
    command_args: dict = None
    file_config: dict = None
    program_path: str = None
    log_path: str = None

    def __init__(self, args: argparse.Namespace = None,
                 base_path: str = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))):
        """
        This method initializes the config manager for VIPER. It saves any arguments given to the program and sets any
        further settings, if necessary.
        Important: command-line settings take precedence over config.json settings!

        :param args: The command-line arguments given to VIPER, as parsed by argparse.ArgumentParser().parse_args()
        :param base_path: The base path of VIPER. Based on this, subfolders will be created and files will be saved.
        """
        if args:
            self.command_args = vars(args)

        try:
            with open(args.config, 'r') as config_file:
                conf_in: dict = json.load(config_file)
                # TODO: Validate file_config? Set defaults? Exit if not all settings are set?

                # recursively flatten JSON structure
                json_conf = {}

                # TODO: Support list entries? Maybe just specify config format to not include them
                def flatten_level(d: dict, add_to: dict, prefix: str = None) -> None:
                    for k, v in d:
                        if isinstance(v, dict):
                            flatten_level(k, add_to, prefix=f"{str(k)}.")
                        else:
                            if prefix:
                                add_to[f"{prefix}{str(k)}"] = v
                            else:
                                add_to[str(k)] = v
                flatten_level(conf_in, json_conf)

                self.file_config = dict(sorted(json_conf))
        except OSError:
            print(f"Could not open config file at '{config_file}'.")
            if "permissive" not in self.command_args:
                print(f"Aborting...")
                sys.exit(1)
            else:
                print(f"Trying to continue with just command line arguments - this might cause errors down the line!")

        self.program_path = base_path

        # Initialize main logger
        if lp := self.get("log_path"):
            self.log_path = self.program_path + os.sep + lp
        else:  # Fallback
            self.log_path = self.program_path + os.sep + "Logs"

        # If the folder where the logs shall be saved doesn't exist, create it
        os.makedirs(self.log_path, exist_ok=True)

        if self.get("verbose"):
            logging.basicConfig(level=logging.DEBUG,
                                filename=os.path.abspath(self.log_path + os.sep + "log.txt"),
                                format="%(asctime)s [%(levelname)s] (%(filename)s:%(module)s): %(message)s",
                                force=True)
        else:
            logging.basicConfig(level=logging.INFO,
                                filename=os.path.abspath(self.log_path + os.sep + "log.txt"),
                                format="%(asctime)s [%(levelname)s] (%(filename)s:%(module)s): %(message)s",
                                force=True)

    def get(self, setting: str) -> Optional[Any]:
        """
        Gets a given setting, if it is set, otherwise returns None. For a list of all setting, see accompanying documentation.
        Important: command-line settings take precedence over config.json settings!

        :param setting: A string representing the name of the setting to be retrieved
        :return: The value of the setting, or None if that setting hasn't been configured
        """
        if self.command_args and setting in self.command_args:
            return self.command_args[setting]
        elif self.file_config and setting in self.file_config:
            return self.file_config[setting]
        else:
            return None

    def log_config(self) -> None:
        """
        Logs all configuration settings to the log file.

        :return: None
        """
        logging.info(f"Running VIPER version {VIPER.VERSION}.")
        logging.info(f"Using the following configuration:")
        if self.file_config:
            logging.info(f"---------- [ CONFIG FILE ] ----------")
            logging.info(pformat(self.file_config))
        if self.command_args:
            logging.info(f"---------- [ COMMAND LINE ARGUMENTS ] ----------")
            logging.info(pformat(self.command_args))
