# VIPER
This is the official repository for the software tool **VIPER: Viral Inhibition via Peptide Engineering and Receptor-Mimicry**.


### Other installations needed
You will need to install [Python 3](https://www.python.org/) (VIPER was developed and tested with Python 3.8), the requirements in the `requirements.txt`, and the [Rosetta Commons](https://www.rosettacommons.org/) software suite.

---
## Quickstart

To reproduce the results from this thesis, follow these steps *after* installing all the requirements outlined above. It's recommended to install the Python packages in a virtual environment.

1. Get the PDB under the accession code **6m0j** and save it at the root level (same folder that `main.py` is in).
2. Open the `config.json` file and add the absolute file path to your RosettaCommons installation under `rosetta_config` → `path`. Also set the amount of cores for Rosetta to use under `rosetta_config` → `use_num_cores`.
3. Run the following command: `python3 -u main.py 6m0j.pdb --vsp_chain E --partner_chain A &`. It is recommended to prefix this with `nohup` if you expect your terminal session to be terminated before VIPER is done and your operating system supports `nohup`. 

You can expect VIPER to be done in less than a day, or even quicker, depending on how many cores you allow it to use.

---
## Running your own workflow

If you simply want to run VIPER with your own PDB file, the setup is very easy. First, make sure that all dependencies are 
installed as outlined above and that the `config.json` is updated correspondingly. Next, identify which chains you want to let VIPER design an inhibitory peptide for.
For the sake of this explanation, let's assume that the viral surface protein has the chain id `V` and the human receptor
protein has the chain id `R`. Save them to `vsp_chain` and `partner_chain`, respectively, in the `config.json`.
Alternatively, you may specify them as command line flags with `--vsp_chain` and `--partner_chain`, respectively.

Next, make sure that all other settings in `config.json` are set to your liking, an explanation of every setting is at
the end of this document.

Finally, you may run VIPER with the command outlined in step 3 in the "Quickstart" section.
This will make VIPER run through the default stages with your PDB file. If you want to adjust this pipeline,
you can adapt it in `VIPER.py`. Currently, the default is as follows:

1. Reduce the PDB to only the relevant chains and renumber the atoms and residues starting from 1.
2. Relax the resultant PDB using Rosetta.
3. Perform a residue energy breakdown of the best relaxed complex using Rosetta.
4. Use this data to derive a candidate peptide (using a selection strategy: `modules` → `stages` → `PeptideGenerator.py` → `SelectionStrategy`).
5. Determine the tertiary structure of this peptide using PEPstrMOD (`modules` → `wrappers` → `PEPstrMODWrapper.py`).
6. Optimize this candidate using the Genetic Algorithm strategy (`modules` → `stages` → `optimize` → `GAStrategy.py`).

If you only want to modify a minor part, VIPER provides some ready-to-use hooks already in the `custom_funcs.py` file. 
There, you can specify a custom function that determines whether a residue should be included when growing subsections
of the partner protein for both the GreedyExpand and FragmentJoiner strategies. You can also implement your own selection
strategy. You may also specify a custom function that gets called when generating mutants during the GA, where you can 
modify the populations as you see fit. Lastly, you can modify how the sSCII value gets combined with the Rosetta score.
If you want to use one of these custom functions, you need to enable the usage in the `config.json`!

If you want to modify / expand VIPER beyond that, you'll have to modify the code directly. `VIPER.py` contains the abstract,
high-level control flow, with most of the logic being implemented in so called 'modules'. There 'wrappers' contain code
to interface with external tools, whereas 'stages' implement internal VIPER logic, such as selecting a subset of residues
from the partner protein. Therefore, your steps should probably be first implementing a wrapper/stage and then tying it into `VIPER.py`.

---
## Support

Are you having issues? Does VIPER crash or not produce the expected output? 
If so, first make sure that all configuration settings are set properly, especially
the PDB chain ids. If that did not fix the issue, you can get more information by adding the `--verbose` command line flag
or setting `"verbose": true` in the config file. Viper should be able to give you more information at the end of the log file
(default: `Logs/log.txt`). If the error message does not help you diagnose the issue or your issue still persists, please
open a GitHub issue, briefly explain the problem you've encountered, and append the last 100 lines of your log file with the `--verbose` flag set.

---

---
## Configuration Options



"log_path": "Logs"
: Which directory to save the logs to. This is also where the short summary report will be saved.

"verbose": false
: Whether to save additional information to the log file. This will increase the size of the log file *significantly*!

"random_seed": 3333333
: Which random seed to use. VIPER will try to use this seed for every tool for which you can set the random seed.

"num_CPU_cores": ""
: If you are using the dynamic concurrent scoring option for the GA, enter how many cores your system has / you want to use.

"results_path": "output"
: Which path to save all generated files besides the logs to.

"vsp_chain": ""
: Which chain in the PDB corresponds to the viral surface protein.

"partner_chain": ""
: Which chain in the PDB corresponds to the partner (receptor) protein.

"reuse_preprocessed": true
: Whether to reuse preprocessed viral surface protein-receptor PDB files, if they are found.

---

"rosetta_config":
: Starts the block of Rosetta config options.

"path": ""
: The path to the RosettaCommons installation. 

"path_out": "rosetta_output"
: Which path to save the output from Rosetta to, unless grouped with other related output. This is mostly where f. e. the
run configurations for Rosetta applications will be saved.
   
"random_seed": -1
: Which random seed to use for Rosetta. If this is set to -1 (default), the general random_seed will be used instead.

"use_num_cores": 8
: How many cores Rosetta should use when running Rosetta applications. 
If this multiplied with your GA population size exceeds the total number of cores available, you *will* starve your system
and overallocate your compute resources!

"archive_intermediate": false
: Whether to archive the intermediate Rosetta output to save space. The code for this exists, but currently isn't integrated.

"prerelax_complex_runs": 100
: How many individual relaxed complexes of the initial PDB to generate.

"relax_partner_runs": 100
: How many individual relaxations of the partner to run. Currently unused.

"relax_peptide_nmr_runs": 40
: How many normal mode relaxations of the peptide to run. Currently unused.

"relax_peptide_bb_runs": 30
: How many backbone relaxations (backrub) of the peptide to run. Currently unused.

"relax_peptide_bb_ntrials": 20000
: How many Monte Carlo trials to run when running backrub relaxation. Currently unused.

"relax_peptide_fast_runs": 30
: How many fast relaxations of the peptide to run. Currently unused.

"docking_runs": 5000
: How many runs of docking to run. Currently unused. Compare to Rosetta `docking_protocol`.

"docking_translation": 3
: How many standard deviations of translation to use when perturbing input structure. Currently unused. Compare to Rosetta `docking_protocol`.

"docking_rotation": 8
: How many standard deviations of rotation to use when perturbing input structure. Currently unused. Compare to Rosetta `docking_protocol`.

"refine_runs": 100
: How many runs of refinement to run on the docked structure(s).

"reb_only_use_best": true
: Whether to only use the best scoring initial relaxed complex for the residue energy breakdown, 
or the average of all relaxations of the initial complex.

---

"gromacs_config"
: Starts the GROMACS config options. Currently unused.

"path": ""
: The path to the GROMACS executable(s). Currently unused.

"random_seed": 0
: What random seed to use. Currently unused.

---

"peptide_generator"
: Starts the config options for the PeptideGenerator.

"use_strategy": "fragment_joiner"
: Which selection strategy to use. Allowed options: `greedy_expand`, `fragment_joiner`.

"linker": "GG"
: What linker to use, if any. Use a string of residues in single letter notation, so 'GG' for two glycines, for example.

"linker_oversize_policy": "truncate"
: Which policy to apply if the linker is bigger than the gap it's supposed to bridge. Options are `truncate` for 
truncating the linker to fit the length by removing amino acids from the end until it matches the gap length, `ignore`
if the linker should always be inserted as is, or `skip` if the linker should be skipped in this case.

"linking_force_length_limit": true
: Whether to also enforce the length limit when linking fragments, i. e. whether linkers should count against the length limit.

"max_length": 18
: The maximum length of the peptide. VIPER will select this many residues from the original partner protein.
If linkers don't count against the length limit, the peptide will likely exceed this number. If linkers do count
against the length limit, the peptide may be shorter than this number.

"reb_energy_cutoff": -0.1
: In absolute terms, what a residue must score better (lower) than in order to be treated favorable and be potentially included in the peptide.

"length_damping": false
: Whether to apply length damping, favoring or penalizing longer fragments.

"length_damping_mode": "quadratic"
: Which damping mode to use. May be either `quadratic` or `linear`.

"length_damping_min_length": 2
: The length after which to apply progressive damping.

"length_damping_max_length": 7
: The length after which to stop applying progressive damping.

"length_damping_initial_mult": 1.0
: The (multiplicative) score modifier to start on at min_length.

"length_damping_final_mult": 0.2
: The (multiplicative) score modifier to end on at max_length.

"length_damping_linear_stepping": ""
: (Optional) The stepping by which to increase / decrease the score per residue length with linear stepping.
Leave empty to let VIPER automatically calculate the stepping that linearily connects the start and end points.

"custom_strategy": false
: Whether to use a custom residue selection strategy (from `custom_funcs.py`).

"greedy_expand": 
: Starts the config options for the GreedyExpand selection strategy.

"max_side_extension": 2
: How far to at most expand a fragment around a start residue in each direction.

"ignore_neighbors": false
: Whether to automatically ignore neighbors and only use the individually highest scoring residues.

"always_include_direct_neighbors": false
: Whether to always include the direct neighbors of a selected residue.

"custom_func": false
: Whether to use a custom inclusion method (from `custom_funcs.py`).

"fragment_joiner":
: Starts the config options for the FragmentJoiner selection strategy.

"use_abs_increase": true
: Whether a residue must increase the total fragment score by a specific absolute amount to be included.

"use_rel_increase": false
: Whether a residue must increase the total fragment score by a specific relative amount to be included.

"mixed_mode_strict": false
: Whether a residue must satisfy both the absolute and relative increase criteria to be included.

"min_abs_increase": -0.2
: The specific absolute value a residue has to improve the score by to meet the abs_increase criterion.

"min_rel_increase": 0.1
: The specific relative percentage a residue has to improve the score by to meet the rel_increase criterion.

"lookahead": 2
: How many residues to look ahead for to see whether a favorable residue can be found once an unfavorable residue has been encountered.

"length_flexibility": 0
: When using the old fragment_combiner method (below), how many residues VIPER is allowed to exceed the length limit by in order to fit in a fragment.

"join_distance_penalty": 4.0
: When using the old fragment_combiner method (below), how far in angstroms fragments may be away from the best fragment, before a score penalty gets applied.

"join_penalty_factor": 0.05
: When using the old fragment_combiner method (below), by how much in percent to reduce the fragment score for each angstrom above the join_distance_penalty it is away.

"fully_join_fragments": false
: Whether to *always* fully join all fragments.

"linker_stretch_factor": 4.0
: How many angstroms another fragment may be away to be elligible for linking for each residue in the linker.

"penalize_lone_residues": true
: Whether to apply a score penalty to fragments consisting of only a single residue.

"lone_residue_penalty": 0.5
: In percent, by how much to reduce the penalty to single-residue fragments.

"pad_lone_residues": true
: Whether to always pad single-residue fragments with neighbors.

"lone_residue_pad_range": 1
: With how many residues in each direction to pad single-residues with.

"old_frag_combiner": false
: Whether to use the old fragment_combiner method, which applies score penalties to fragments based on their distance to the best scoring fragment.

"custom_func": false
: Whether to use a custom inclusion method (from `custom_funcs.py`).

---

"pepstrmod_config":
: Starts the config options for the PEPstrMODWrapper.

"email": " "
: Which email address to fill into the email address field of the job submission form. If it is to be left empty, put in a space!

"simulation_time": "100ps"
: For how long to run the MD simulation for. May be either "100ps" or "50ps".

"environment": "vac"
: In which environment to simulate the peptide. 

"download_topology": true
: Whether to make the MD sim topology files available on the results page. These will *not* be automatically downloaded by VIPER!

"cluster_analysis": false
: Whether to perform a cluster analysis.

"download_trajectory": true
: Whether to make the MD sim trajectory files available on the results page. This will *not* be automatically downloaded by VIPER!

"do_energy_rms_graph": true
: Whether to show an energy RMS graph on the results page. This will *not* be automatically downloaded by VIPER!

"wait_interval": 60
: For how many seconds to wait between checking if the results are ready.

---

"optimize"
: Starts the config options for the optimization stage.

"ga"
: Starts the config options for the genetic algorithm strategy.

"select_percent": 0.3
: How many parents to select from the previous generation to generate offspring for the next generation, in percent.

"selection_mode": "ROULETTEWHEEL"
: How to select parents. Available options are "ROULETTEWHEEL" for a weighted random selection based on their fitness, 
"BESTONLY" to only select the x% best individuals as parents, or "UNIFORM" to randomly select parents from the previous generation.

"selection_with_replacement": true
: Whether to sample parents from the previous generation with replacement. This may cause a parent to appear multiple 
times in the next generation.

"crossover_mode": "MULTIPLE"
: How often to cross over between parents when generating offspring. "SINGLE" only has a single crossover point, whereas
"MULTIPLE" allows for multiple crossover points to appear.

"crossover_chance": 0.1
: The chance, in percent, at each residue to cross over between parents.

"mutation_rate": 0.05
: The chance, in percent, at each residue to mutate into another residue.

"mutation_bias": "BLOSUM62_shifted"
: How to bias the residue mutation. May be the name of any of the matrices present in the folder `util` → `substitution_matrices` folder.
Use `UNIFORM` for unbiased mutation.

"num_generations": 5
: How many generations to run the genetic algorithm for.

"join_pops_after": -1
: After how many generations to merga all populations into one. Set to -1 to never merge populations.

"getstruc_backoff": 600
: Range, in seconds, to select a random wait time for before submitting the individual (peptide sequence) to the tertiary structure prediction service.
This staggers the submissions in order to not overload / stress the service too much, if it is an external one, such as PEPstrMOD.

"num_relax_individual": 10
: How many relaxations to run of the predicted tertiary structure of the peptide.

"dynamic_concurrent_scoring": false
: Whether to use a dynamic concurrent scoring paradigm, where VIPER automatically tries to determine the best number of cores to use
for each run of Rosetta during scoring, so as not to overallocate cores.

"contact_checking"
: Starts the config options for contact checking during scoring.

"adjust_score": false
: Whether to use the results of contact checking to adjust the score of the individual.

"emit_warning": true
: Whether to save the results of the contact checking to disk.

"mismatch_tolerance": 2
: How many original contacts may be lost per residue before applying a score penalty.

"nearby_partner_tolerance": 1
: How many residues in each direction a contact may still count as an original one. I. e., the original contact was 
residue 42 on chain X and this setting equals 2, a new contact with either residue 40, 41, 42, 43, or 44 still counts as
an original contact.

"penalty": 0.02
: In percent, by how much to reduce the individual's score by for each residue that misses enough original contacts.

"scii"
: Starts the config options for sSCII calculation during scoring.

"adjust_score": true
: Whether to adjust the individual's score based on the sSCII value.

"radius": 7
: The radius, in angstrom, to use for sSCII calculation.

"threshold": 0.4063
: The threshold to distinguish stable/unstable residues by for sSCII calculation.

"stepping_width": 0.1
: For how long to remain on a modifier value before moving on to the next modifier value.

"bonus": 0.05
: By how much to increase / decrease the modifier value per step.

"custom_func": false
: Whether to use a custom sSCII bases score modifier method (from `custom_funcs.py`).

"custom_addin_mutate": false
: Whether to apply a custom mutation function in addition to the existing one (from `custom_funcs.py`)

