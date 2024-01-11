from __future__ import annotations

import copy
import itertools
import logging
import math
import os
import pprint
import random
from pathlib import Path
from typing import List, Callable, Union

import ConfigManager
from . import OptimizationStrategy
from modules.wrappers import RosettaWrapper
from modules.wrappers.PEPstrMODWrapper import PEPstrMODWrapper
from util import SCII, file_utils, PDBtool
from util.BLOSUM import BLOSUM

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

    def select(self, population: Population, n: int = None):
        take_num = math.ceil(self.config["select_percent"] * len(population)) if n is None else n
        lookup = {}
        # Get scores for population
        for individual in population:
            if individual in self.score_repo:
                lookup[individual] = self.score_repo[individual]
                continue
            sc = self._score_func(individual)
            lookup[individual] = sc
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
                    if gene1 is None:
                        flip = not flip
                    else:
                        offspring.append(gene1)
                elif not flip:
                    offspring.append(gene2)
        if self.config["crossover_mode"] == "SINGLE":
            crossover_point = random.randrange(len(combined))
            offspring += parent1[:crossover_point]
            offspring += parent2[crossover_point:]
        return offspring

    def mutate(self, individual):
        _ = copy.deepcopy(individual)
        for n, gene in enumerate(individual):
            if random.random() < self.config["mutation_rate"]:
                if bmatrix := self.config["mutation_bias"]:
                    _[n] = random.choices(population=list(bmatrix[gene].keys()), weights=list(bmatrix[gene].values()),
                                          k=1)
                else:
                    _[n] = random.choices(population=list(BLOSUM.BLOSUM62_shifted.get(gene)), k=1)
        return _

    def score(self, peptide: str) -> float:
        # Get 3D structure
        pdb = PEPstrMODWrapper.submit_peptide(sequence=peptide)
        use_path = os.path.join(self.out_path, f"gen{self.generation}_{peptide}")
        pdb = file_utils.make_file(os.path.normpath(use_path + "base.pdb"), pdb)

        # Superimpose onto receptor
        pdb, rms = PDBtool.superimpose_single(pdb, self.ref,
                                              query_chain=f"{PDBtool.get_chains(os.path.normpath(pdb))[0]}",
                                              ref_chain=f"{cm().get('partner_chain')}",
                                              out_path=[use_path, "aligned", f"{peptide}_aligned.pdb"])

        # Get binding energy
        self.rw.run(RosettaWrapper.Flags().interface_analyzer, options={
            "-in:file:s": pdb,
            "-out:file:score_only": os.path.join(use_path, "interface_score.sc"),
        })
        best_pdb, scores = RosettaWrapper.ScoreFileParser.get_extremum(
            os.path.normpath(
                os.path.join(use_path, "interface_score.sc")),
            "dG_separated")
        best_pdb = os.path.normpath(use_path + best_pdb + ".pdb")
        best_rosetta_score = scores["dG_separated"]

        # Calculate SCII
        scii = SCII.scii_for_pdb(best_pdb, radius=self.config["score_scii_radius"])

        # Calculate composite score
        # Since 0.3 - 0.4 was the area of ambiguity stated in the original SCII paper, with larger values indicating
        # more stable peptides and smaller values indicating less stable peptides
        # Therefore a percentage based bonus or malus will be applied if the value is higher or lower thant 0.3 - 0.4
        bonus = round((scii - 0.35), 1) * 0.5 + 1
        self.score_repo[peptide] = best_rosetta_score * bonus
        return self.score_repo[peptide]

    def run(self):
        for i in range(self.config["num_generations"]):
            logging.debug(f"Genetic Algorithm, generation {self.generation}")
            for pop in self.populations:
                logging.debug(f"Old pop: {pprint.pformat(pop)}")
                old_num = len(pop)
                families = []
                parents = self._select(pop)
                # Since we need to pair up parents, check for even number of parents. If the number of parents is odd,
                # duplicate the first entry and append it.
                if len(parents) % 2 != 0:
                    parents.append(parents[0])
                parents = [(parents[i], parents[i + 1]) for i in range(0, len(parents), 2)]
                for p1, p2 in parents:
                    families.append(p1)
                    families.append(p2)
                    offspring = self._crossover(p1, p2)
                    if offspring not in families:
                        families.append(offspring)
                muts = []
                for member in families:
                    muts.append(self._mutate(member))
                # Pad to previous length with mutants of randomly chosen finalists of this generation
                while len(muts) < old_num:
                    candidate = self._mutate(random.choice(muts))
                    if candidate not in muts:
                        muts.append(candidate)
                pop.update(muts)
                logging.debug(f"New pop: {pop}")
                curr_best = pop.get_asc()[0]
                logging.info(
                    f"Generation {self.generation}, best candidate: {curr_best} ({self.score_repo[curr_best]})")
