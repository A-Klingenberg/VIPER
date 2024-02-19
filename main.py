"""
Main file for VIPER.
"""
import argparse
import os

from ConfigManager import ConfigManager
from VIPER import VIPER


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    # TODO: Add all config file arguments here as well
    parser.add_argument("PDB", help="The path to a PDB file containing the viral surface protein bound to the human "
                                    "cell surface protein", type=str)
    parser.add_argument("--config", help="The path to a config file", type=str)
    parser.add_argument("--log_path", help="The path to where log files should be written", type=str)
    parser.add_argument("--results_path", help="The path to where the output of all files should be written",
                        type=str)
    parser.add_argument("--verbose", help="Set this flag to log additional information", action="store_true")
    parser.add_argument("--permissive", help="Set this flag to have VIPER continue running even if it encounters "
                                             "problems which might lead to unexpected behaviour", type=bool,
                        default=True)
    parser.add_argument("--num_CPU_cores", help="How many CPU cores to use", type=int)
    parser.add_argument("--vsp_chain", help="The chain id of the viral surface protein", type=str)
    parser.add_argument("--partner_chain", help="The chain id of the partner protein, for which the"
                                                "inhibitory peptide candidates shall be generated", type=str)
    # rosetta config
    parser.add_argument("--rosetta_config.path", help="The path to your rosetta executable", type=str)
    parser.add_argument("--rosetta_config.path_out",
                        help="The path to where the output of rosetta runs shall be written", type=str)
    parser.add_argument("--rosetta_config.relax_protein_runs", help="How many protein relaxation runs should be run",
                        type=int)
    parser.add_argument("--rosetta_config.relax_xml_runs", help="How many xml relaxation runs should be run", type=int)
    parser.add_argument("--rosetta_config.relax_bb_runs", help="How many backbone relaxation runs should be run",
                        type=int)
    parser.add_argument("--rosetta_config.relax_fast_runs", help="How many fast relaxation runs should be run",
                        type=int)
    parser.add_argument("--rosetta_config.docking_runs", help="How many docking runs should be run", type=int)
    parser.add_argument("--rosetta_config.refine_runs", help="How many refinement runs should be run", type=int)


    # gromacs config

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    # noinspection PyArgumentList
    ConfigManager(force_refresh_singleton=True, args=args, base_path=os.path.abspath(os.path.dirname(os.path.realpath(__file__))))
    v = VIPER(ConfigManager.get_instance().get("PDB"))
    v.run()


if __name__ == "__main__":
    main()
