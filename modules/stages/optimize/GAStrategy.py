from __future__ import annotations

import copy
import itertools
import logging
import math
import multiprocessing
import os
import pprint
import random
import shutil
import time
from pathlib import Path
from typing import List, Callable, Union, Tuple

import ConfigManager
from modules.wrappers import RosettaWrapper
from modules.wrappers.PEPstrMODWrapper import PEPstrMODWrapper
from util import SCII, file_utils, PDBtool
from util.BLOSUM import BLOSUM
from . import OptimizationStrategy

cm = ConfigManager.ConfigManager.get_instance


# Genetic Algorithm Optimization Strategy

class Population:
    individuals: List[str] = []
    score_func: Callable = None

    def __init__(self, score_func: Callable, individuals: List[str] = None):
        self.score_func = score_func
        if individuals is not None:
            self.individuals = individuals

    def get_desc(self) -> List:
        return sorted(self.individuals, key=self.score_func if self.score_func is not None else None, reverse=True)

    def get_asc(self) -> List:
        return sorted(self.individuals, key=self.score_func if self.score_func is not None else None, reverse=False)

    def update(self, individuals: List) -> Population:
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
        return f"[{', '.join(self.individuals)}"[:-2] + "]"

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
    populations: List[Population] = None
    _score_func: Callable = None
    _select: Callable = None
    _crossover: Callable = None
    _mutate: Callable = None
    config: dict = None
    score_repo: dict = None
    generation: int = 0
    rw: RosettaWrapper = None
    out_path: Path = None
    ref: Union[Path, str] = None
    vsp: Union[Path, str] = None

    def __init__(self, ref_pdb: Union[Path, str], populations: List[Population], score_func: Callable = None,
                 select: Callable = None,
                 crossover: Callable = None, mutate: Callable = None, config: dict = None,
                 propagate_score_func: bool = True, metric: str = "MIN"):
        super().__init__()
        self.ref = ref_pdb
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
        self.config["mutation_bias"] = config.get("mutation_bias", BLOSUM.BLOSUM62_shifted)
        self.config["num_generations"] = config.get("num_generations", 5)
        self.config["getstruc_backoff"] = config.get("getstruc_backoff",
                                                     4 * 60)  # wait 0-4 minutes, try to not stress webservice
        self.score_repo = {}
        self.generation = 0
        self.rw = RosettaWrapper.RosettaWrapper()
        self.out_path = Path(os.path.join(cm().get("results_path"), "GA"))
        if metric not in ["MIN", "MAX"]:
            self.metric = "MIN"
        else:
            self.metric = metric
        if propagate_score_func:
            for pop in self.populations:
                pop.score_func = self._score_func
        os.makedirs(os.path.join(cm().get("results_path"), "GA"), exist_ok=True)
        self.vsp = PDBtool.remove_chain(self.ref, [cm().get("partner_chain")],
                                       os.path.join(cm().get("results_path"), "GA", "vsp.pdb"))

    def select(self, population: Population, n: int = None):
        take_num = math.ceil(self.config["select_percent"] * len(population)) if n is None else n
        # Get scores for population
        # Scoring may be very time intensive, so do this concurrently
        with multiprocessing.Pool() as pool:
            pool.map(self._score_func, population)
        # We can now be sure that we have scores for every individual in this population
        lookup = {individual: self.score_repo[individual]["total"] for individual in population}
        if self.config["selection_mode"] == "UNIFORM":
            return random.choices(list(lookup.keys()), k=take_num)
        # Order by fitness
        ordered = None
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
            return random.choices([ind for ind, _ in ordered], weights=[fit for _, fit in ordered], k=take_num)
        else:
            # Return whole population
            return [ind for ind, _ in ordered]

    def crossover(self, parent1, parent2):
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
        _ = list(copy.deepcopy(individual))
        for n, gene in enumerate(individual):
            if random.random() < self.config["mutation_rate"]:
                if bmatrix := self.config["mutation_bias"]:
                    _[n] = random.choices(population=list(bmatrix[gene].keys()), weights=list(bmatrix[gene].values()),
                                          k=1)[0]
                else:
                    _[n] = random.choices(population=list(BLOSUM.BLOSUM62_shifted.get(gene)), k=1)[0]
        return "".join(_)

    def _getstruc(self, peptide: str) -> Tuple[Path, Path]:
        if backoff := self.config["getstruc_backoff"]:
            time.sleep(random.randrange(0, backoff))
        pdb = PEPstrMODWrapper.submit_peptide(sequence=peptide)
        use_path_items = ["GA", f"gen{self.generation}_{peptide}"]
        pdb = file_utils.make_file(["GA", f"gen{self.generation}_{peptide}", "base.pdb"], pdb)

        # Superimpose onto receptor and create merged pdb
        peptide_pdb, rms = PDBtool.superimpose_single(pdb, self.ref,
                                                      query_chain=f"{PDBtool.get_chains(os.path.normpath(pdb))[0]}",
                                                      ref_chain=f"{cm().get('partner_chain')}",
                                                      out_path=[*use_path_items, f"{peptide}_aligned.pdb"])
        return PDBtool.join(peptide_pdb, self.vsp), peptide_pdb

    def score(self, peptide: str) -> float:
        if peptide in self.score_repo:
            return self.score_repo[peptide]["total"]
        start = time.time()

        # Get 3D structure
        complex_pdb, peptide_pdb = self._getstruc(peptide)

        # Make sure complex for binding energy scoring is energetically favorable / relaxed
        relax_path = os.path.join(complex_pdb.parent, "relax")
        os.makedirs(os.path.join(relax_path, "complex"), exist_ok=True)
        self.rw.run(RosettaWrapper.Flags().relax_base, options={
            "-in:file:s": complex_pdb,
            "-nstruct": 100,
            "-out:path:all": os.path.join(relax_path, "complex"),
            "-out:suffix": "_relax",
        })
        best_complex, complex_scores = RosettaWrapper.ScoreFileParser.get_extremum(
            os.path.normpath(os.path.join(relax_path, "complex", f"score_relax.sc")), "total_score")
        shutil.copyfile(os.path.join(relax_path, "complex", best_complex + ".pdb"),
                        os.path.join(complex_pdb.parent, "best_complex.pdb"))

        # Get binding energy
        score_path = os.path.join(self.out_path, f"gen{self.generation}_{peptide}", "interface_score.sc")
        self.rw.run(RosettaWrapper.Flags().interface_analyzer, options={
            "-in:file:s": False,
            "-s": os.path.join(complex_pdb.parent, "best_complex.pdb"),
            "-out:file:score_only": os.path.normpath(score_path),
        })
        best_pdb, scores = RosettaWrapper.ScoreFileParser.get_extremum(
            os.path.normpath(score_path),
            "dG_separated")
        best_pdb = f"GA/gen{self.generation}_{peptide}/{best_pdb}.pdb"
        best_rosetta_score = scores["dG_separated"]

        # Calculate SCII
        scii = SCII.scii_for_pdb(peptide_pdb, radius=self.config["score_scii_radius"])

        # Calculate composite score
        # Since 0.3 - 0.4 was the area of ambiguity stated in the original SCII paper, with larger values indicating
        # more stable peptides and smaller values indicating less stable peptides
        # Therefore a percentage based bonus or malus will be applied if the value is higher or lower thant 0.3 - 0.4
        bonus = round((scii - 0.35), 1) * 0.5 + 1
        self.score_repo[peptide] = {"total": best_rosetta_score * bonus, "rosetta_score": best_rosetta_score,
                                    "SCII": scii, "calc_scii_bonus": bonus}
        logging.debug(f"Scoring {peptide} took {(time.time() - start) / 60:.1f} minutes!")
        return self.score_repo[peptide]["total"]

    def run(self):
        for i in range(self.config["num_generations"]):
            logging.debug(f"Genetic Algorithm, generation {self.generation}")
            pop_count = 0
            for pop in self.populations:
                logging.debug(f"Old pop [#{pop_count}] : {pprint.pformat(pop)}")
                old_num = len(pop)
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
                pop.update(finalists)
                logging.debug(f"New pop [#{pop_count}] : {pop}")
                curr_best = pop.get_asc()[0]
                logging.info(
                    f"Generation {self.generation}, best candidate: {curr_best} ({self.score_repo[curr_best]['total']})")
                pop_count += 1
            self.generation += 1
        best = min(self.score_repo.items(), key=lambda item: item[1]["total"])
        logging.info(f"Best found candidate is {best[0]} with score {best[1]['total']}.")
        logging.debug(pprint.pformat(self.score_repo))
