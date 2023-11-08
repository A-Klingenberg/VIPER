import os.path
from pathlib import Path
from typing import Union

import ConfigManager
from modules.wrappers import RosettaWrapper
from util import file_utils

cm = ConfigManager.ConfigManager.get_cm


# TODO: Does this need to be a class?
class VIPER:
    clean_pdb_reference = None
    predocked_reference_pdb = None
    rw = None

    def __init__(self, pdb: Union[str, Path]):
        self.rw = RosettaWrapper.RosettaWrapper()
        self.predocked_reference_pdb = Path(pdb)
        pass

    def run(self):
        """
        This method shall execute each computational step needed to generate an inhibitory peptide.

        :return:
        """


        pass

    def preprocess_pdb(self):

        pass

    def do_energy_breakdown(self):
        # prerelax to identify significant residues

        out_path = self.rw.make_dir(["energy_breakdown", "prerelax"])
        self.rw.run(RosettaWrapper.Flags.relax_complex_for_REB, options={
            "-in:file:s": os.path.normpath(self.predocked_reference_pdb),
            "-out:path:all": out_path,
        })
        relaxed_initials = file_utils.gather_files(out_path)
        out_path = os.path.normpath(os.path.join(out_path, ".."))
        for n, pdb in enumerate(relaxed_initials, start=1):
            self.rw.run(RosettaWrapper.Flags.residue_energy_breakdown, options={
                "-in:file:s": pdb,
                "-out:file:silent": os.path.join(out_path, f"energy_breakdown_{Path(pdb).name}_{n}.out"),
            })
        energy_breakdowns = file_utils.gather_files(out_path, type="out")
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
