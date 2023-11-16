import logging
import os.path
import shutil
from pathlib import Path
from typing import Union, List

import ConfigManager
from modules.stages import PeptideGenerator
from modules.wrappers import RosettaWrapper
from modules.wrappers.PEPstrMODWrapper import PEPstrMODWrapper
from util import file_utils, PDBtool

cm = ConfigManager.ConfigManager.get_cm


# TODO: Does this need to be a class?
class VIPER:
    base_pdb = None
    reference_renum_pdb = None
    reference_renum_relaxed = None
    rw = None
    selection_strat = None
    candidate_counter = 0

    def __init__(self, pdb: Union[str, Path]):
        self.rw = RosettaWrapper.RosettaWrapper()
        self.base_pdb = Path(pdb)
        self.reference_renum_relaxed = Path("6m0j_renum_relaxed.pdb")
        self.selection_strat = PeptideGenerator._SelectionStrategies.get_strategy()()
        self.candidate_counter = 0
        pass

    def run(self):
        """
        This method shall execute each computational step needed to generate an inhibitory peptide.

        :return:
        """
        self.preprocess_pdb()
        nlist = self.do_energy_breakdown()
        candidate = self.generate_peptide(nlist)
        peptide_structure = self.get_tertiary_structure(candidate,
                                                        ["candidates", str(self.candidate_counter), f"candidate.pdb"])
        p, _ = PDBtool.superimpose_single(os.path.normpath(peptide_structure),
                                          os.path.normpath(self.reference_renum_relaxed),
                                          query_chain=f"{PDBtool.get_chains(os.path.normpath(peptide_structure))[0]}",
                                          ref_chain=f"{cm().get('partner_chain')}",
                                          out_path=["candidates", str(self.candidate_counter), f"superimposed.pdb"])
        curr_candidate_dir = Path(os.path.join(cm().get("results_path"), "candidates", str(self.candidate_counter)))
        # Do a fast relaxation of the peptide while it is pinned in place
        self.rw.run(RosettaWrapper.Flags.relax_complex_for_REB.copy(), options={
            "-in:file:s": os.path.normpath(p),
            "-out:path:all": os.path.join(curr_candidate_dir, "relax"),
            "-in:file:native": os.path.normpath(self.reference_renum_pdb),
        })
        # Do residue energy breakdown to get binding energy
        self.rw.run(RosettaWrapper.Flags.residue_energy_breakdown.copy(), options={
            "-in:file:l": os.path.normpath(
                file_utils.make_pdb_ensemble_list(os.path.join(curr_candidate_dir, "relax"),
                                                  os.path.join(curr_candidate_dir,
                                                               f"candidate_{self.candidate_counter}_relax_ensemble"))),
            "-out:file:silent": os.path.join(curr_candidate_dir,
                                             f"energy_breakdown_candidate{self.candidate_counter}_ensemble.out")
        })

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
            "-nstruct": cm().get("rosetta_config.prerelax_complex_runs"),
        })
        initials = file_utils.make_pdb_ensemble_list(intermediary_dir,
                                                     os.path.normpath(
                                                         os.path.join(intermediary_dir, "..",
                                                                      "reference_ensemble")))
        logging.debug(f"Relaxed base PDB and found files '{initials}'")
        best_pdb, scores = RosettaWrapper.ScoreFileParser.get_extremum(
            os.path.normpath(
                os.path.join(intermediary_dir,
                             f"score{RosettaWrapper.Flags.relax_complex_for_REB['-out:suffix']}.sc")),
            "total_score")
        logging.debug(f"Best PDB is '{best_pdb}' with score {scores['total_score']}")
        self.reference_renum_relaxed = Path(os.path.normpath(
            os.path.join(intermediary_dir, "..", self.base_pdb.name[:-4] + "_renum_relaxed.pdb")))
        shutil.copyfile(os.path.join(intermediary_dir, best_pdb + ".pdb"),
                        os.path.join(intermediary_dir, self.reference_renum_relaxed))
        logging.info(f"Finished preparing PDB! Prepared PDB is saved to '{self.reference_renum_relaxed}'")

    def do_energy_breakdown(self):
        # Compare residue involvement in all the relaxations of the experimental structure
        out_path = os.path.normpath(self.rw.make_dir(["residue_energy_breakdown", "reference_pdb"]))
        self.rw.run(RosettaWrapper.Flags.residue_energy_breakdown.copy(), options={
            "-in:file:l": os.path.normpath(
                os.path.join(cm().get("results_path"), "reference", "reference_ensemble")),
            "-out:file:silent": os.path.join(out_path,
                                             f"energy_breakdown_{Path(self.reference_renum_relaxed).name[:-4]}.out")
        })
        return RosettaWrapper.REBprocessor.process_multipose(
            os.path.join(out_path, f"energy_breakdown_{Path(self.reference_renum_relaxed).name[:-4]}.out"))

    def generate_peptide(self, nlist: List[RosettaWrapper.REBprocessor.Node]) -> List[RosettaWrapper.REBprocessor.Node]:
        candidate = self.selection_strat.reduce(from_chain=cm().get("partner_chain"),
                                                to_chain=cm().get("vsp_chain"),
                                                nodes=nlist)
        return candidate

    def get_tertiary_structure(self, sequence: List[RosettaWrapper.REBprocessor.Node],
                               save_out: List[Union[str, Path]]) -> Path:
        peptide_sequence = "".join([PDBtool.three_to_one(n.amino_acid) for n in sequence])
        structure_pdb = PEPstrMODWrapper.submit_peptide(sequence=peptide_sequence)
        return file_utils.make_file(path=save_out, content=structure_pdb)

    def prepare_docking(self):
        pass

    def dock_peptide(self):
        pass

    def restore_pdb_formatting(self):
        pass

    def calc_stability(self):
        # MD?
        pass
