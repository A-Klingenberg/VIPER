import logging
import os.path
import shutil
from pathlib import Path
from typing import Union

import ConfigManager
from modules.wrappers import RosettaWrapper
from util import file_utils, PDBtool

cm = ConfigManager.ConfigManager.get_cm


# TODO: Does this need to be a class?
class VIPER:
    base_pdb = None
    reference_renum_pdb = None
    reference_renum_relaxed = None
    rw = None
    chains = {}

    def __init__(self, pdb: Union[str, Path]):
        self.rw = RosettaWrapper.RosettaWrapper()
        self.base_pdb = Path(pdb)
        pass

    def run(self):
        """
        This method shall execute each computational step needed to generate an inhibitory peptide.

        :return:
        """
        self.preprocess_pdb()
        self.do_energy_breakdown()

        pass

    def preprocess_pdb(self):
        logging.info(f"Preprocessing PDB '{self.base_pdb}'...")
        intermediary_dir = Path(os.path.join(cm().get("results_path"), "reference", "intermediary"))
        os.makedirs(intermediary_dir, exist_ok=True)
        self.reference_renum_pdb = PDBtool.renumber_ascending(os.path.normpath(self.base_pdb),
                                                              os.path.normpath(os.path.join(intermediary_dir,
                                                                                            "..",
                                                                                            self.base_pdb.name[:-4],
                                                                                            "_renum.pdb")))
        logging.debug(f"Relaxing renumbered PDB...")
        self.rw.run(RosettaWrapper.Flags.relax_complex_for_REB, options={
            "-in:file:s": os.path.normpath(self.reference_renum_pdb),
            "-out:path:all": intermediary_dir,
            "-in:file:native": os.path.normpath(self.base_pdb),
        })
        initials = file_utils.make_pdb_ensemble_list(intermediary_dir,
                                                     os.path.normpath(os.path.join(intermediary_dir, "..")))
        logging.debug(f"Relaxed base PDB and found files '{initials}'")
        best_pdb, scores = RosettaWrapper.ScoreFileParser.get_extremum(
            os.path.normpath(os.path.join(intermediary_dir, "score.sc")),
            "total_score")
        logging.debug(f"Best PDB is '{best_pdb}' with score {scores['total_score']}")
        self.reference_renum_relaxed = Path(
            os.path.join(intermediary_dir, self.base_pdb.name[:-4] + "_relaxed_reference.pdb"))
        shutil.copyfile(os.path.join(intermediary_dir, best_pdb + ".pdb"),
                        os.path.join(intermediary_dir, self.reference_renum_relaxed))
        logging.info(f"Finished preparing PDB! Prepared PDB is saved to '{self.reference_renum_relaxed}'")

    def do_energy_breakdown(self):
        # Compare residue involvement in all the relaxations of the experimental structure
        relaxed_initials = file_utils.gather_files(os.path.join(cm().get("results_path"), "reference", "intermediary"))
        out_path = os.path.normpath(self.rw.make_dir(["residue_energy_breakdown", "base_pdb"]))
        for n, pdb in enumerate(relaxed_initials, start=1):
            self.rw.run(RosettaWrapper.Flags.residue_energy_breakdown, options={
                "-in:file:s": os.path.normpath(pdb),
                "-out:file:silent": os.path.join(out_path, f"energy_breakdown_{Path(pdb).name}_{n}.out"),
            })
        energy_breakdowns = file_utils.gather_files(out_path, filetype="out")
        nlists = []
        for breakdown in energy_breakdowns:
            nlists.append(RosettaWrapper.REBprocessor.read_in(breakdown)[0])
        primary = []
        for n, nlist in enumerate(nlists):
            if n == 0:
                primary = nlist
            else:
                for n_index, node in enumerate(nlist):
                    for chain, energy in node.strength.items():
                        if chain in primary[n_index].strength.items():
                            primary[n_index].strength[chain] += energy
                        else:
                            primary[n_index].strength[chain] = energy
        for node in primary:
            for chain, energy in node.strength.items():
                node.strength[chain] /= len(nlists)
        return primary, nlists[0]

    def generate_peptide(self):
        pass

    def get_tertiary_structure(self):
        pass

    def prepare_docking(self):
        pass

    def dock_peptide(self):
        pass

    def restore_pdb_formatting(self):
        pass

    def calc_stability(self):
        # MD?
        pass
