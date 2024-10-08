# VIPER   [![DOI](https://zenodo.org/badge/682657296.svg)](https://zenodo.org/doi/10.5281/zenodo.10897858)
This is the official repository for the software tool **VIPER: Viral Inhibition via Peptide Engineering and Receptor-Mimicry**.


### Other installations needed
You will need to install [Python 3](https://www.python.org/) (VIPER was developed and tested with Python 3.8), the requirements in the `requirements.txt`, and the [Rosetta Commons](https://www.rosettacommons.org/) software suite.

---
## Manual

VIPER is extensively documented. You can have a look at the Quickstart, descriptions of the internal algorithms, or a glossary of all config options in the [VIPER manual](VIPER_Manual.pdf). It is highly recommended to review the manual first.

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
