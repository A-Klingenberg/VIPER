import copy
import logging
import os.path
import pprint
import random
import shutil
import sys
from pathlib import Path
from typing import Union, List, Any, Tuple

import ConfigManager
from modules.stages import PeptideGenerator
from modules.stages.optimize.GAStrategy import GAStrategy, Population
from modules.wrappers import RosettaWrapper
from modules.wrappers.PEPstrMODWrapper import PEPstrMODWrapper
from util import file_utils, PDBtool
from util.substitution_matrices import submat

cm = ConfigManager.ConfigManager.get_instance


# TODO: Does this need to be a class?
class VIPER:
    base_pdb = None
    base_pdb_cleaned = None
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
        self.summary_logger = logging.getLogger("summary")
        self.summary_logger.info(f"Trying to generate peptides for {self.base_pdb.name}!")

    def run(self) -> None:
        """
        Runs all steps necessary to generate an inhibitory peptide.

        :return: None
        """

        self.preprocess_pdb()
        reb_nodes = self.do_energy_breakdown()
        raw_candidate = self.generate_peptide(reb_nodes)
        curr_candidate_dir = Path(os.path.join(cm().get("results_path"), "candidates", str(self.candidate_counter)))
        peptide_structure = self.get_tertiary_structure(raw_candidate,
                                                        [curr_candidate_dir, f"candidate.pdb"])
        best, score = self.optimize(raw_candidate, pad_to=10, out_path=curr_candidate_dir)

    def preprocess_pdb(self) -> Path:
        """
        Preprocesses the PDB file by renumbering atom and residues ids from 1 and relaxing the complex using Rosetta.

        :return: The path to the relaxed and renumbered PDB file
        """
        logging.info(f"Preprocessing PDB '{self.base_pdb}'...")
        intermediary_dir = Path(os.path.join(cm().get("results_path"), "reference", "intermediary"))
        reuse = os.path.abspath(os.path.join(intermediary_dir, "..", self.base_pdb.name[:-4] + "_renum_relaxed.pdb"))
        if os.path.isfile(reuse):
            logging.warning(f"Already found a preprocessed file, potentially from an earlier run: {reuse}")
            if cm().get("reuse_preprocessed"):
                logging.warning("Reusing that file!")
                self.summary_logger.info(f"Reusing a preprocessed reference PDB ({reuse})!")
                self.reference_renum_relaxed = Path(os.path.normpath(
                    os.path.join(intermediary_dir, "..", self.base_pdb.name[:-4] + "_renum_relaxed.pdb")))
                cm().ref_relax = self.reference_renum_relaxed
                return self.reference_renum_relaxed
            else:
                logging.warning("Ignoring and potentially overwriting that file! ...")
        os.makedirs(intermediary_dir, exist_ok=True)
        chains = PDBtool.get_chains(self.base_pdb)
        rem_chains = []
        for c in chains:
            if c not in cm().get("vsp_chain") and c not in cm().get("partner_chain"):
                rem_chains.append(c)
        self.base_pdb_cleaned = PDBtool.remove_chain(self.base_pdb, rem_chains, os.path.join(intermediary_dir, "..",
                                                                                             self.base_pdb.name[
                                                                                             :-4] + "_cleaned.pdb"))
        self.base_pdb_cleaned = PDBtool.reorder_chains(self.base_pdb_cleaned,
                                                       chain_order=f"{cm().get('partner_chain')}{cm().get('vsp_chain')}")
        self.reference_renum_pdb = PDBtool.renumber_ascending(os.path.normpath(self.base_pdb_cleaned),
                                                              os.path.normpath(os.path.join(intermediary_dir,
                                                                                            "..",
                                                                                            self.base_pdb_cleaned.name[
                                                                                            :-4] +
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
        cm().ref_relax = self.reference_renum_relaxed
        shutil.copyfile(os.path.join(intermediary_dir, best_pdb + ".pdb"),
                        os.path.join(intermediary_dir, self.reference_renum_relaxed))
        logging.info(f"Finished preparing PDB! Prepared PDB is saved to '{self.reference_renum_relaxed}'")
        self.summary_logger.info(f"Input PDB ({self.base_pdb.name}) was successfully preprocessed! The reference PDB "
                                 f"that will be used going forward is saved under '{self.reference_renum_relaxed}'.")
        return self.reference_renum_relaxed

    def do_energy_breakdown(self) -> List[RosettaWrapper.REBprocessor.Node]:
        """
        Runs the Rosetta Energy Breakdown and parses the output. May read in a multipose output file, for more info
        refer to RosettaWrapper.REBprocessor.process_multipose

        :return: A list of RosettaWrapper.REBprocessor.Node objects holding the energy info for each residue
        """
        # Compare residue involvement in all the relaxations of the experimental structure
        out_path = os.path.normpath(self.rw.make_dir(["residue_energy_breakdown", "reference_pdb"]))
        if cm().get("rosetta_config.reb_only_use_best"):
            self.rw.run(RosettaWrapper.Flags().residue_energy_breakdown, options={
                "-in:file:s": os.path.normpath(self.reference_renum_relaxed),
                "-out:file:silent": os.path.join(out_path,
                                                 f"energy_breakdown_{Path(self.reference_renum_relaxed).name[:-4]}.out")
            })
        else:
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
        self.summary_logger.info(f"An initial candidate peptide has been generated: "
                                 f"{' '.join([_.__repr__() for _ in candidate])[:-1]}")
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
        result = file_utils.make_file(path=save_out, content=structure_pdb)
        curr_candidate_dir = Path(os.path.join(cm().get("results_path"), "candidates", str(self.candidate_counter)))
        # Do a fast relaxation of the peptide while it is pinned in place
        self.rw.run(RosettaWrapper.Flags().relax_pinned_positions, options={
            "-in:file:s": os.path.normpath(result),
            "-out:path:all": os.path.join(curr_candidate_dir, "relax"),
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
        self.summary_logger.info(f"Successfully got tertiary structure for sequence {peptide_sequence}! "
                                 f"It is saved to '{result}'.")
        return result

    def prepare_docking(self):
        pass

    def dock_peptide(self):
        pass

    def restore_pdb_formatting(self):
        pass

    def calc_stability(self):
        # MD?
        pass

    def optimize(self, raw_candidate: List[RosettaWrapper.REBprocessor.Node], initpops: List[Population] = None,
                 pad_to: int = -1, out_path: Union[str, Path] = None) -> Tuple[Any, float]:
        """
        Optimizes the passed raw candidate. You can change this out for whatever optimization procedure you like.

        :param raw_candidate: A list of REBprocessor Nodes that represents the candidate you want to optimize.
        :param initpops: If you already have a set of populations, you may pass it here. If this isn't set, you must set
            pad_to to an integer > 1.
        :param pad_to: If you don't have an initial population to optimize, set this number to configure the number of
            random mutants to generate based on the raw candidate to form the initial population. If this isn't set, you
            must pass a list of populations in initpops.
        :param out_path: An optional path to a folder that the GA results should be saved to.
        :return: A tuple of the best candidate and its aggregate score.
        """
        if initpops is None and pad_to <= 1:
            if not cm().get("permissive", default=True):
                raise ValueError("You must pass either an initial population or configure the number of individuals to "
                                 "expand the raw candidate to!")
            else:
                pad_to = 9
        bmatrix = submat.SubMat(cm().get("optimize.ga.mutation_bias", "BLOSUM62_shifted.json"))

        if initpops is not None:
            ga = GAStrategy(ref_pdb=self.reference_renum_relaxed, populations=initpops,
                            orig_pep_contacts=raw_candidate, out_base=out_path,
                            config={
                                "mutation_bias": bmatrix,
                            })
        else:
            def mutate(individual):
                _ = list(copy.deepcopy(individual))
                for n, gene in enumerate(_):
                    if random.random() < cm().get("optimize.ga.mutation_rate", 0.05):
                        _[n] = random.choices(population=list(bmatrix.get(gene).keys()),
                                              weights=list(bmatrix.get(gene).values()),
                                              k=1)[0]
                return "".join(_)

            pepseq = "".join([PDBtool.three_to_one(n.amino_acid) for n in raw_candidate])

            initpop = []
            while len(initpop) < pad_to - 1:  # One less, because we need to append raw_candidate as well
                mut = mutate(pepseq)
                if mut not in initpop and mut != pepseq:
                    initpop.append(mut)
            initpop.append(pepseq)
            ga = GAStrategy(ref_pdb=self.reference_renum_relaxed, populations=[Population(initpop)],
                            orig_pep_contacts=raw_candidate, out_base=out_path,
                            config={
                                "mutation_bias": bmatrix,
                            })

        best_candidate, best_score = ga.run()
        return best_candidate, best_score
