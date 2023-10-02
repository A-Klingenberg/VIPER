import ConfigManager

VERSION = 0.1
configmanager: ConfigManager = None


# TODO: Does this need to be a class?
class VIPER:
    clean_pdb_reference = None
    predocked_reference_pdb = None

    def __init__(self, cm: ConfigManager):
        self.cm = cm
        global configmanager
        configmanager = cm

    def preprocess_pdb(self):
        pass

    def do_energy_breakdown(self):
        pass

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

    def run(self):
        """
        This method shall execute each computational step needed to generate an inhibitory peptide.

        :return:
        """
        pass
