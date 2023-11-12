import copy
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
                                                                                            self.base_pdb.name[:-4] +
                                                                                            "_renum.pdb")))
        logging.debug(f"Relaxing renumbered PDB...")
        self.rw.run(RosettaWrapper.Flags.relax_complex_for_REB.copy(), options={
            "-in:file:s": os.path.normpath(self.reference_renum_pdb),
            "-out:path:all": intermediary_dir,
            "-in:file:native": os.path.normpath(self.reference_renum_pdb),
        })
        initials = file_utils.make_pdb_ensemble_list(intermediary_dir,
                                                     os.path.normpath(
                                                         os.path.join(intermediary_dir, "..", "reference_ensemble")))
        logging.debug(f"Relaxed base PDB and found files '{initials}'")
        best_pdb, scores = RosettaWrapper.ScoreFileParser.get_extremum(
            os.path.normpath(
                os.path.join(intermediary_dir, f"score{RosettaWrapper.Flags.relax_complex_for_REB['-out:suffix']}.sc")),
            "total_score")
        logging.debug(f"Best PDB is '{best_pdb}' with score {scores['total_score']}")
        self.reference_renum_relaxed = Path(os.path.normpath(
            os.path.join(intermediary_dir, "..", self.base_pdb.name[:-4] + "_renum_relaxed.pdb")))
        shutil.copyfile(os.path.join(intermediary_dir, best_pdb + ".pdb"),
                        os.path.join(intermediary_dir, self.reference_renum_relaxed))
        logging.info(f"Finished preparing PDB! Prepared PDB is saved to '{self.reference_renum_relaxed}'")

    def do_energy_breakdown(self):
        # Compare residue involvement in all the relaxations of the experimental structure
        relaxed_initials = file_utils.gather_files(os.path.join(cm().get("results_path"), "reference", "intermediary"))
        out_path = os.path.normpath(self.rw.make_dir(["residue_energy_breakdown", "reference_pdb"]))
        self.rw.run(RosettaWrapper.Flags.residue_energy_breakdown.copy(), options={
            "-in:file:l": os.path.normpath(os.path.join(cm().get("results_path"), "reference", "reference_ensemble")),
            "-out:file:silent": os.path.join(out_path,
                                             f"energy_breakdown_{Path(self.reference_renum_relaxed).name[:-4]}.out")
        })

        poses = RosettaWrapper.REBprocessor.read_in(
            os.path.join(out_path, f"energy_breakdown_{Path(self.reference_renum_relaxed).name[:-4]}.out"))
        aggregate = []
        backup = []
        first = True
        for pose, nlist in poses.items():
            if first:
                # Add dummy records
                aggregate = [RosettaWrapper.REBprocessor.Node(n.amino_acid, n.residue_id, n.chain, n.partners.copy(),
                                                              copy.deepcopy(n.strength)) for n in nlist]
                backup = [RosettaWrapper.REBprocessor.Node(n.amino_acid, n.residue_id, n.chain, n.partners.copy(),
                                                           copy.deepcopy(n.strength)) for n in nlist]
                first = False
            else:
                for n_index, node in enumerate(nlist):
                    for chain, energy in node.strength.items():
                        if chain in aggregate[n_index].strength.items():
                            aggregate[n_index].strength[chain] += energy
                        else:
                            aggregate[n_index].strength[chain] = energy
        use_pose = poses.popitem(last=False)
        for n_index, node in enumerate(aggregate):
            for chain, energy in node.strength.items():
                node.strength[chain] = node.strength[chain] / len(poses)

        return aggregate, backup

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
