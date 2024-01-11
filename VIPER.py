import copy
import logging
import os.path
import random
import shutil
import sys
from pathlib import Path
from typing import Union, List

import Bio

import ConfigManager
from modules.stages import PeptideGenerator
from modules.stages.optimize.GAStrategy import GAStrategy, Population
from modules.wrappers import RosettaWrapper
from modules.wrappers.PEPstrMODWrapper import PEPstrMODWrapper
from util import file_utils, PDBtool
from util.BLOSUM import BLOSUM

cm = ConfigManager.ConfigManager.get_instance


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
        self.selection_strat = PeptideGenerator._SelectionStrategies.get_strategy()()
        self.candidate_counter = 0
        pass

    def run(self) -> None:
        """
        Runs all steps necessary to generate an inhibitory peptide.

        :return: None
        """
        self.preprocess_pdb()
        nlist = self.do_energy_breakdown()
        candidate = self.generate_peptide(nlist)

        # Generate _without_ linkers!
        peptide_structure = self.get_tertiary_structure([n for n in candidate if n.orig_res_id is not None],
                                                        ["candidates", str(self.candidate_counter),
                                                         f"candidate.pdb"])
        # peptide_structure = self.get_tertiary_structure(candidate,
        #                                                ["candidates", str(self.candidate_counter), f"candidate.pdb"])

        curr_candidate_dir = Path(os.path.join(cm().get("results_path"), "candidates", str(self.candidate_counter)))

        def mutate(individual):
            _ = copy.deepcopy(individual)
            for n, gene in enumerate(individual):
                if random.random() < 0.05:
                    if bmatrix := self.config["mutation_bias"]:
                        _[n] = random.choices(population=list(bmatrix[gene].keys()),
                                              weights=list(bmatrix[gene].values()),
                                              k=1)
                    else:
                        _[n] = random.choices(population=list(BLOSUM.BLOSUM62.keys()), k=1)
            return _

        pepnodes = [n for n in candidate if n.orig_res_id is not None]
        pepseq = "".join([PDBtool.three_to_one(n.amino_acid) for n in pepnodes])

        initpop = []
        while len(initpop) < 10:
            mut = mutate(pepseq)
            if mut not in initpop:
                initpop.append(mut)
        initpop.append(pepseq)
        ga = GAStrategy(self.reference_renum_pdb, [Population(lambda _: _, initpop)])
        ga.run()


        """
        # Do a fast relaxation of the peptide while it is pinned in place
        self.rw.run(RosettaWrapper.Flags().relax_pinned_positions, options={
            "-in:file:s": os.path.normpath(peptide_structure),
            "-out:path:all": os.path.join(curr_candidate_dir, "relax"),
            "-in:file:native": os.path.normpath(self.reference_renum_pdb),
        })
        # Do residue energy breakdown to get binding energy
        self.rw.run(RosettaWrapper.Flags().residue_energy_breakdown, options={
            "-in:file:l": os.path.normpath(
                file_utils.make_pdb_ensemble_list(os.path.join(curr_candidate_dir, "relax"),
                                                  os.path.join(curr_candidate_dir,
                                                               f"candidate_{self.candidate_counter}_relax_ensemble"))),
            "-out:file:silent": os.path.join(curr_candidate_dir,
                                             f"energy_breakdown_candidate{self.candidate_counter}_ensemble.out")
        })
        """

    def preprocess_pdb(self) -> Path:
        """
        Preprocesses the PDB file by renumbering atom and residues ids from 1 and relaxing the complex using Rosetta.

        :return: The path to the relaxed and renumbered PDB file
        """
        logging.info(f"Preprocessing PDB '{self.base_pdb}'...")
        intermediary_dir = Path(os.path.join(cm().get("results_path"), "reference", "intermediary"))
        os.makedirs(intermediary_dir, exist_ok=True)
        self.reference_renum_pdb = PDBtool.renumber_ascending(os.path.normpath(self.base_pdb),
                                                              os.path.normpath(os.path.join(intermediary_dir,
                                                                                            "..",
                                                                                            self.base_pdb.name[:-4] +
                                                                                            "_renum.pdb")))
        logging.debug(f"Relaxing renumbered PDB...")
        self.rw.run(RosettaWrapper.Flags().relax_pinned_positions, options={
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
                             f"score{RosettaWrapper.Flags().relax_pinned_positions['-out:suffix']}.sc")),
            "total_score")
        logging.debug(f"Best PDB is '{best_pdb}' with score {scores['total_score']}")
        self.reference_renum_relaxed = Path(os.path.normpath(
            os.path.join(intermediary_dir, "..", self.base_pdb.name[:-4] + "_renum_relaxed.pdb")))
        shutil.copyfile(os.path.join(intermediary_dir, best_pdb + ".pdb"),
                        os.path.join(intermediary_dir, self.reference_renum_relaxed))
        logging.info(f"Finished preparing PDB! Prepared PDB is saved to '{self.reference_renum_relaxed}'")
        return self.reference_renum_relaxed

    def do_energy_breakdown(self) -> List[RosettaWrapper.REBprocessor.Node]:
        """
        Runs the Rosetta Energy Breakdown and parses the output. May read in a multipose output file, for more info
        refer to RosettaWrapper.REBprocessor.process_multipose

        :return: A list of RosettaWrapper.REBprocessor.Node objects holding the energy info for each residue
        """
        # Compare residue involvement in all the relaxations of the experimental structure
        out_path = os.path.normpath(self.rw.make_dir(["residue_energy_breakdown", "reference_pdb"]))
        self.rw.run(RosettaWrapper.Flags().residue_energy_breakdown, options={
            "-in:file:l": os.path.normpath(
                os.path.join(cm().get("results_path"), "reference", "reference_ensemble")),
            "-out:file:silent": os.path.join(out_path,
                                             f"energy_breakdown_{Path(self.reference_renum_relaxed).name[:-4]}.out")
        })
        return RosettaWrapper.REBprocessor.process_multipose(
            os.path.join(out_path, f"energy_breakdown_{Path(self.reference_renum_relaxed).name[:-4]}.out"))

    def generate_peptide(self, nlist: List[RosettaWrapper.REBprocessor.Node]) -> List[RosettaWrapper.REBprocessor.Node]:
        """
        Generates a peptide based on the residue energy breakdown output. For more info refer to the selection
        strategies.

        :param nlist: A list of RosettaWrapper.REBprocessor.Node objects, the result of reading in an energy breakdown
            file of a complex
        :return: A list of RosettaWrapper.REBprocessor.Node objects representing the generated peptide
        """
        candidate = self.selection_strat.reduce(from_chain=cm().get("partner_chain"),
                                                to_chain=cm().get("vsp_chain"),
                                                nodes=nlist)
        return candidate

    def get_tertiary_structure(self, sequence: List[RosettaWrapper.REBprocessor.Node],
                               save_out: List[Union[str, Path]]) -> Path:
        """
        Gets the tertiary structure for a passed sequence of residues (Nodes)

        :param sequence: A list of RosettaWrapper.REBprocessor.Node objects representing the sequence of the peptide
        :param save_out: Where to save the PDB with the tertiary structure to.
        :return: A path object to the saved PDB with the tertiary structure.
        """
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
