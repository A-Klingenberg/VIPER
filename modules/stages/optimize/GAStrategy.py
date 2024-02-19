from __future__ import annotations

import copy
import itertools
import json
import logging
import math
import multiprocessing
import operator
import os
import pprint
import random
import shutil
import time
from pathlib import Path
from typing import List, Callable, Union, Tuple, Any

import ConfigManager
from modules.wrappers import RosettaWrapper
from modules.wrappers.PEPstrMODWrapper import PEPstrMODWrapper
from util import SCII, file_utils, PDBtool
from util.substitution_matrices import submat
from . import OptimizationStrategy

cm = ConfigManager.ConfigManager.get_instance


class Population:
    """
    A population hold multiple individuals. Individuals within a population are meant to compete against each other and
    possibly reproduce with other individuals within this population during regular generations of the genetic algorithm

    This implementation of a population assumes that the individual('s genome) can be represented as a str!
    """
    individuals: List[str] = []
    score_func: Callable = None
    score_repo: dict = None

    def __init__(self, individuals: List, score_func: Callable = None, scores: dict = None):
        self.individuals = individuals
        self.score_func = score_func
        self.score_repo = scores

    def get_desc(self, scores: dict = None) -> List:
        """
        Gets the individuals in descending order of their scores.

        :param scores: (Optional) A dict of scores to use. If this is not None, every individual MUST have a numerical
            entry in the scores dictionary!
        :return: The ordered list of individuals, from high scores to low scores.
        """
        use_scores = scores if scores else self.score_repo
        return sorted(self.individuals, key=lambda i: self.score_func(i, shared_dict=use_scores, verbosity=cm().get(
            "verbose"))["total"] if self.score_func is not None else None, reverse=True)

    def get_asc(self, scores=None) -> List:
        """
        Gets the individuals in ascending order of their scores.

        :param scores: (Optional) A dict of scores to use. If this is not None, every individual MUST have a numerical
            entry in the scores dictionary!
        :return: The ordered list of individuals, from low scores to high scores.
        """
        use_scores = scores if scores else self.score_repo
        return sorted(self.individuals, key=lambda i: self.score_func(i, shared_dict=use_scores, verbosity=cm().get(
            "verbose"))["total"] if self.score_func is not None else None, reverse=False)

    def update(self, individuals: List) -> Population:
        """
        Replaces the individuals of this population.

        :param individuals: A new list of individuals to replace the current population with.
        :return: The population itself
        """
        self.individuals = individuals
        return self

    def __len__(self):
        return len(self.individuals)

    def __getitem__(self, item):
        if type(item) is int or type(item) is slice:
            return self.individuals[item]
        else:
            raise TypeError("Can only index into individuals using integers or slices!")

    def __iter__(self):
        return (_ for _ in self.individuals)

    def __str__(self):
        """
        The representation of the population is a list of the str representation of its individuals

        :return: A str representation of the population
        """
        return f"[{', '.join(self.individuals)}" + "]"

    def __repr__(self):
        return self.__str__()

    # TODO: This is probably too hacky, write something more robust
    def __contains__(self, item):
        # Assume individual is a string of amino acids in single-letter abbreviation
        for ind in self.individuals:
            if str(item) == str(ind):
                return True
        return False


class GAStrategy(OptimizationStrategy.OptimizationStrategy):
    """
    This class manages runs of a genetic algorithm optimization strategy.
    Has support for multiple concurrent populations.
    """
    populations: List[Population] = None
    _score_func: Callable = None
    _select: Callable = None
    _crossover: Callable = None
    _mutate: Callable = None
    config: dict = None
    score_repo: dict = None
    generation: int = 0
    rw: RosettaWrapper = None
    out_path: Union[Path, str] = None
    ref: Union[Path, str] = None
    vsp: Union[Path, str] = None
    _cust_addin_mutate: bool = False
    do_contacts_check: bool = False

    def __init__(self, ref_pdb: Union[Path, str], populations: List[Population], ref_reb: Union[Path, str] = None,
                 orig_pep_contacts: List[RosettaWrapper.REBprocessor.Node] = None, score_func: Callable = None,
                 select: Callable = None, crossover: Callable = None, mutate: Callable = None, config: dict = None,
                 propagate_score_func: bool = True, metric: str = "MIN"):
        """
        Initializes the genetic algorithm strategy with the passed parameters.

        :param ref_pdb: A path to the reference pdb for which to optimize the peptide.
        :param populations: A list of populations to optimize. If you only want to use a single Population, you still
            have to wrap it in a list! ( [...] )
        :param ref_reb: A reference residue energy breakdown of the original VSP+receptor complex to use in scoring.
            (Optional) If this argument isn't given, original contact checking will not be performed.
        :param orig_pep_contacts: The starting peptide to use for residue-residue contact analysis. (Optional) If this
            argument isn't given, original contact checking will not be performed.
        :param score_func: (Optional) The score function to use on the individual. If not supplied, uses the builtin
        :param select: (Optional) The function to select parents. If not supplied, uses the builtin
        :param crossover: (Optional) The function to apply the crossover operation. If not supplied, uses the builtin
        :param mutate: (Optional) The function to apply the mutation operation. If not supplied, uses the builtin
        :param config: (Optional) Any settings / parameters of the genetic algorithm you want to change from the default
        :param propagate_score_func: Whether to propagate the score function supplied here to the populations
            (default: True)
        :param metric: Whether a smaller ("MIN") or larger ("MAX") score is better (default: "MIN")
        """
        super().__init__()
        self.ref = ref_pdb
        self.ref_reb = RosettaWrapper.REBprocessor.process_multipose(ref_reb)
        self.orig_pep_contacts = orig_pep_contacts
        self.contact_dict = {}
        if self.ref_reb and self.orig_pep_contacts:
            self.do_contacts_check = True
            for n, node in enumerate(self.orig_pep_contacts):
                if not isinstance(node, RosettaWrapper.REBprocessor.Node):
                    if cm().get("permissive"):
                        logging.warning(f"The orig_pep_contacts that was passed doesn't exclusively contain "
                                        f"RosettaWrapper.REBprocessor.Node objects! Skipping entry {n} ({node}).")
                    else:
                        raise ValueError("The orig_pep_contacts that was passed doesn't exclusively contain "
                                         "RosettaWrapper.REBprocessor.Node objects!")
                # Only save interactions with nodes that aren't on chain and preferable for binding (<0)
                self.contact_dict[n] = [partner[0] for partner in node.partners if partner[0][-1] != node.chain
                                        and partner[1] < 0.0] if node.orig_res_id is not None else []
        self._score_func = score_func if score_func is not None else self.score
        self.populations = populations
        self._select = select if select is not None else self.select
        self._crossover = crossover if crossover is not None else self.crossover
        self._mutate = mutate if mutate is not None else self.mutate
        self.config = config if config is not None else {}
        self.config["select_percent"] = config.get("select_percent", 0.3)
        self.config["selection_mode"] = config.get("selection_mode", "ROULETTEWHEEL")
        self.config["crossover_mode"] = config.get("crossover_mode", "MULTIPLE")
        self.config["crossover_chance"] = config.get("crossover_chance", 0.1)
        self.config["mutation_rate"] = config.get("mutation_rate", 0.05)
        self.config["mutation_bias"] = config.get("mutation_bias", submat.UNIFORM)
        self.config["num_generations"] = config.get("num_generations", 5)
        self.config["getstruc_backoff"] = config.get("getstruc_backoff",
                                                     10 * 60)  # wait 0-10 minutes, try to not stress webservice
        self.config["num_relax_individual"] = config.get("num_relax_individual", 10)
        self.score_repo = {}
        self.generation = 0
        self.rw = RosettaWrapper.RosettaWrapper()
        self.out_path = os.path.normpath(os.path.join(cm().get("results_path"), "GA"))
        if metric not in ["MIN", "MAX"]:
            self.metric = "MIN"
        else:
            self.metric = metric
        if propagate_score_func:
            for pop in self.populations:
                pop.score_func = self._score_func
        os.makedirs(self.out_path, exist_ok=True)
        self.vsp = PDBtool.remove_chain(self.ref, [cm().get("partner_chain")],
                                        os.path.join(cm().get("results_path"), "GA", "vsp.pdb"))
        self._cust_addin_mutate = cm().get("optimize.ga.custom_addin_mutate")
        logging.info(f"Instantiating with following parameters: {pprint.pformat(self.config, compact=True)}")

    def select(self, population: Population, n: int = None) -> List:
        """
        Implements the selection operator. Chooses parents from the passed population and returns them in a list.

        :param population: The population to select parents for.
        :param n: (Optional) How many parents to select. If this is not specified, it derives the number from the config
        :return: A list of parent individuals
        """
        take_num = round(self.config["select_percent"] * len(population)) if n is None else n
        if take_num == 0:
            take_num = 1
        # restrict view to current population
        lookup = {individual: self.score_repo.get(individual, self._score_func(individual))["total"] for individual in
                  population}
        if self.config["selection_mode"] == "UNIFORM":
            return random.choices(list(lookup.keys()), k=take_num)
        # Order by fitness
        if self.metric == "MIN":
            ordered = sorted(lookup.items(), key=lambda tup: tup[1])
        elif self.metric == "MAX":
            ordered = sorted(lookup.items(), key=lambda tup: tup[1], reverse=True)
        else:
            # standard: sort with lower fitness being best
            ordered = sorted(lookup.items(), key=lambda tup: tup[1])
        # Select individuals
        if self.config["selection_mode"] == "BESTONLY":
            return [ind for ind, _ in ordered[:take_num]]
        if self.config["selection_mode"] == "ROULETTEWHEEL":
            return random.choices([ind for ind, _ in ordered], weights=[abs(fit) for _, fit in ordered], k=take_num)
        else:
            # Return whole population
            return [ind for ind, _ in ordered]

    def crossover(self, parent1, parent2):
        """
        Apply the crossover operation to the two parents and returns the new individual. Assumes the genome is a str of
        amino acids in single letter notation.

        :param parent1: The first parent
        :param parent2: The second parent
        :return: A new individual based on the parents
        """
        parent1, parent2 = (parent1, parent2) if len(parent1) < len(parent2) else (parent2, parent1)
        combined = itertools.zip_longest(parent1, parent2, fillvalue=None)
        offspring = []
        if self.config["crossover_mode"] == "MULTIPLE":
            flip = True
            end = False
            for gene1, gene2 in combined:
                if random.random() < self.config["crossover_chance"] and not end:
                    flip = not flip
                if flip:
                    if gene1 is None:  # handle if one individual has a shorter genome than the other
                        flip = not flip
                    else:
                        offspring.append(gene1)
                elif not flip:
                    offspring.append(gene2)
        if self.config["crossover_mode"] == "SINGLE":
            crossover_point = random.randrange(max(len(parent1), len(parent2)))
            offspring += parent1[:crossover_point]
            offspring += parent2[crossover_point:]
        return offspring

    def mutate(self, individual):
        """
        Applies the mutation operator. Can be biased in its choice of mutation based on BLOSUM matrices. Assumes the
        genome is a string amino acids in single letter notation. Potentially returns the same individual, if the
        mutation rate is fairly low.

        :param individual: The individual to mutate
        :return: The (possibly) mutated individual.
        """
        _ = list(copy.deepcopy(individual))
        for n, gene in enumerate(individual):
            if random.random() < self.config["mutation_rate"]:
                _[n] = random.choices(population=list(self.config["mutation_bias"].get(gene).keys()),
                                      weights=list(self.config["mutation_bias"].get(gene).values()),
                                      k=1)[0]
        return "".join(_)

    def _getstruc(self, peptide: str) -> Tuple[Path, Path]:
        """
        Helper function that takes in a peptide as a str of amino acids in single letter notation and returns a path to
        a PDB of the complex of VSP + superimposed peptide and a path to a PDB of the predicted tertiary structure of
        the PDB.

        :param peptide: The peptide as a str of amino acids in single letter notation
        :return: A tuple of (Path to PDB of VSP+superimposed peptide complex, Path to PDB of peptide tertiary structure)
        """
        if backoff := self.config["getstruc_backoff"]:
            # Wait a bit to stagger requests to webservice
            delay = random.randrange(0, backoff)
            logging.debug(f"Delaying {delay} seconds before getting structure for {peptide}")
            time.sleep(delay)
        pdb = PEPstrMODWrapper.submit_peptide(sequence=str(peptide))
        use_path_items = ["GA", f"gen{self.generation}_{peptide}"]
        pdb = file_utils.make_file(["GA", f"gen{self.generation}_{peptide}", "base.pdb"], pdb)

        # Superimpose onto receptor and create merged pdb
        peptide_pdb, rms = PDBtool.superimpose_single(pdb, self.ref,
                                                      query_chain=f"{PDBtool.get_chains(os.path.normpath(pdb))[0]}",
                                                      ref_chain=f"{cm().get('partner_chain')}",
                                                      out_path=[*use_path_items, f"{peptide}_aligned.pdb"])
        return PDBtool.join(peptide_pdb, self.vsp), peptide_pdb

    def __getstate__(self):
        temp_dict = self.__dict__.copy()
        del temp_dict["rw"]  # Remove RosettaWrapper from pickled GAStrategy object to hopefully stop RecursionErrors
        del temp_dict["ref_reb"]  # same reasoning here
        return temp_dict

    def score(self, peptide: str, shared_dict=None, base_log_path: Union[Path, str] = None,
              verbosity: bool = False) -> Tuple[Any, dict]:
        """
        Scores a peptide and updates the shared_dict (or self.score_repo if not specified). Logs to a separate file
        in the base_log_path given the verbosity. Returns a dict of the scores, with the value under "total" being the
        aggregate score.

        :param peptide: The peptide as a str of amino acids in single letter notation
        :param shared_dict: (Optional) A dictionary shared between runs of this method. If not specified,
            defaults to self.score_repo - You should either make sure this is a DictProxy for safe concurrent sharing,
            or a JSON string that can be deserialized to a dictionary. You can NOT assume that changes made to the
            dictionary will be transparent/permanent outside of this method! This is intended as one-way, read-only
            sharing of values.
        :param base_log_path: The base log path to put the log file under. Will be written to "score_log.txt" in the
            subfolder "gen{self.generation}_{peptide}"
        :param verbosity: Whether to enable verbose logging. (default: False)
        :return: A tuple of an individual and it's score. The score should be a dictionary with the individual terms,
            and a total score under the key 'total'
        """
        if shared_dict is None:
            shared_dict = self.score_repo
        elif not isinstance(shared_dict, dict) and isinstance(shared_dict, str):
            shared_dict = json.loads(shared_dict)
        if base_log_path is None:
            base_log_path = self.out_path
        if peptide in shared_dict:
            # This case should never happen, since this is already checked before calling score. This is a safety
            return shared_dict[peptide]["total"]
        # Because we are likely running this in a separate process, log to separate files for a more robust approach
        os.makedirs(os.path.normpath(os.path.join(base_log_path, f"gen{self.generation}_{peptide}")), exist_ok=True)
        if verbosity:
            logging.basicConfig(level=logging.DEBUG,
                                filename=os.path.normpath(
                                    os.path.join(base_log_path, f"gen{self.generation}_{peptide}", "score_log.txt")),
                                format="%(asctime)s [%(levelname)s] (%(filename)s:%(module)s): %(message)s",
                                force=True)
        else:
            logging.basicConfig(level=logging.INFO,
                                filename=os.path.normpath(
                                    os.path.join(base_log_path, f"gen{self.generation}_{peptide}", "score_log.txt")),
                                format="%(asctime)s [%(levelname)s] (%(filename)s:%(module)s): %(message)s",
                                force=True)
        start = time.time()

        # Get 3D structure
        complex_pdb, peptide_pdb = self._getstruc(peptide)

        # Make sure complex for binding energy scoring is energetically favorable / relaxed
        relax_path = os.path.join(complex_pdb.parent, "relax")
        os.makedirs(os.path.join(relax_path, "complex"), exist_ok=True)
        RosettaWrapper.RosettaWrapper().run(RosettaWrapper.Flags().relax_base, flag_suffix=peptide, options={
            "-in:file:s": complex_pdb,
            "-nstruct": self.config["num_relax_individual"],
            "-out:path:all": os.path.join(relax_path, "complex"),
            "-out:suffix": "_relax",
        })
        best_complex, complex_scores = RosettaWrapper.ScoreFileParser.get_extremum(
            os.path.normpath(os.path.join(relax_path, "complex", f"score_relax.sc")), "total_score")
        shutil.copyfile(os.path.join(relax_path, "complex", best_complex + ".pdb"),
                        os.path.join(complex_pdb.parent, "best_complex.pdb"))

        # Get interface energy
        score_path = os.path.join(self.out_path, f"gen{self.generation}_{peptide}", "interface_score.sc")
        RosettaWrapper.RosettaWrapper().run(RosettaWrapper.Flags().interface_analyzer, flag_suffix=peptide, options={
            "-in:file:s": False,
            "-s": os.path.join(complex_pdb.parent, "best_complex.pdb"),
            "-out:file:score_only": os.path.normpath(score_path),
        })
        best_pdb, scores = RosettaWrapper.ScoreFileParser.get_extremum(
            os.path.normpath(score_path),
            "dG_separated")
        best_rosetta_score = scores["dG_separated"]

        score_modifications = []
        if len(self.contact_dict) != 0:
            logging.debug("Doing score modification based on original contacts / interactions")
            # Get binding energy
            score_path = os.path.join(self.out_path, f"gen{self.generation}_{peptide}", "reb_score.sc")
            RosettaWrapper.RosettaWrapper().run(RosettaWrapper.Flags().residue_energy_breakdown,
                                                     flag_suffix=peptide, options={
                "-in:file:s": os.path.join(complex_pdb.parent, "best_complex.pdb"),
                "-out:file:silent": score_path,
            })
            interactions = RosettaWrapper.REBprocessor.process_multipose(score_path)
            mismatch_tolerance = 0
            # Every inter-residue interaction that doesn't appear for the corresponding residue in the new candidate
            # that exceeds the tolerance (# of mismatches) incurs a 5% score penalty
            for i in range(len(self.orig_pep_contacts)):
                new_partners = [partner[0] for partner in interactions[i].partners]
                ref_partners = self.contact_dict[i]
                mismatch_count = 0
                logging.debug(f"Doing residue {i} with partners {pprint.pformat(new_partners)}. "
                              f"Original: {pprint.pformat(ref_partners)}")
                for partner in ref_partners:
                    if partner not in new_partners:
                        mismatch_count += 1
                        if mismatch_count > mismatch_tolerance:
                            score_modifications.append((operator.mul, 0.95))
                            break

        # Calculate SCII
        scii = SCII.scii_for_pdb(peptide_pdb, radius=self.config["score_scii_radius"])

        # Calculate composite score
        # Since 0.3 - 0.4 was the area of ambiguity stated in the original SCII paper, with larger values indicating
        # more stable peptides and smaller values indicating less stable peptides
        # Therefore a percentage based bonus or malus will be applied if the value is higher or lower thant 0.3 - 0.4
        bonus = round((scii - 0.35), 1) * 0.5 + 1
        score_modifications.append((operator.mul, bonus))

        final_score = best_rosetta_score
        for mod in score_modifications:
            final_score = mod[0](final_score, mod[1])

        tup = (peptide, {"total": best_rosetta_score * bonus, "rosetta_score": best_rosetta_score,
                         "SCII": scii, "calc_scii_bonus": bonus})

        logging.debug(
            f"Scoring {peptide} (score {tup[1]['total']}) took {(time.time() - start) / 60:.1f} minutes!")
        return tup

    def run(self) -> Tuple[Any, float]:
        """
        Runs the genetic algorithm for optimizing the peptide. Returns a tuple of the best individual as a str of amino
         acids in single letter noration and the associated aggregate score.

        :return: A tuple of (best individual, aggregate score)
        """
        for i in range(self.config["num_generations"]):
            logging.info(f"Genetic Algorithm, generation {self.generation}")
            pop_count = 0
            for pop in self.populations:
                logging.info(f"Old pop [#{pop_count}] : {pprint.pformat(pop)}")
                old_num = len(pop)
                # Get scores for population
                # Size task chunks correctly, i. e. don't oversubscribe processors
                # First, determine how many cores to use per tasklet from config, or use the
                # next highest power of two of 20% of the number of processors in the system
                # Second, see how often that number fits into the total number of processors to be used (default: on the
                # system and use as many workers (or +1 if necessary) to ensure that the system isn't
                # oversubscribed to cores
                num, extra = divmod(cm().get("rosetta_config.use_num_cores",
                                             2**(math.ceil(0.2 * os.cpu_count())-1).bit_length()),
                                    cm().get("num_CPU_cores", os.cpu_count()))
                if extra != 0:
                    num += 1
                logging.debug(f"Determined chunksize ({num}) with ({cm().get('num_CPU_cores', os.cpu_count())}) processors")
                # Scoring may be very time intensive, so do this concurrently for every individual within this
                # population for which we don't already have a score
                with multiprocessing.Pool(cm().get("num_CPU_cores", os.cpu_count())) as pool:
                    print(f"{[(individual, '{}', self.out_path, cm().get('verbose')) for individual in pop if individual not in self.score_repo]}")
                    for result in pool.starmap(self._score_func,
                                               [(individual, {}, self.out_path, cm().get("verbose")) for
                                                individual in pop if individual not in self.score_repo], chunksize=num):
                        self.score_repo[result[0]] = result[1]
                # We can now be sure that we have scores for every individual in this population
                families = []
                parents = self._select(pop)
                # Since we need to pair up parents, check for even number of parents. If the number of parents is odd,
                # duplicate the first entry and append it.
                if len(parents) % 2 != 0:
                    parents.append(parents[0])
                parents_paired = [(parents[i], parents[i + 1]) for i in range(0, len(parents), 2)]
                for p1, p2 in parents_paired:
                    offspring = self._crossover(p1, p2)
                    while offspring == p1 or offspring == p2:
                        offspring = self._crossover(p1, p2)
                    if offspring not in families:
                        families.append(offspring)
                finalists = []
                for member in families:
                    finalists.append(self._mutate(member))
                # Pad to previous length with mutants of randomly chosen offspring of this generation
                while len(finalists) < old_num - len(parents):
                    candidate = self._mutate(random.choice(families))
                    if candidate not in finalists:
                        finalists.append(candidate)
                for p in parents:
                    finalists.append(p)
                curr_best = min({i: self.score_repo[i]['total'] for i in pop}.items(), key=lambda t: t[1])
                logging.info(
                    f"Generation {self.generation}, best candidate: {curr_best[0]} ({curr_best[1]})")
                pop.update(finalists)
                if self._cust_addin_mutate:
                    import custom_funcs
                    add_to = custom_funcs.addin(self, pop)
                    for individual in add_to:
                        pop.individuals.append(individual)
                logging.info(f"New pop [#{pop_count}] : {pop}")
                pop_count += 1
            self.generation += 1
        best = min(self.score_repo.items(), key=lambda item: item[1]["total"])
        logging.info(f"Best found candidate is {best[0]} with score {best[1]['total']}.")
        logging.debug(f"Dumping score repo: {pprint.pformat(self.score_repo)}")
        return best[0], best[1]["total"]
