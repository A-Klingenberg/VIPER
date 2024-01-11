# This code has been adapted from previous work by Austin Seamann, Dario Ghersi, and Ryan Ehrlich.
# Please refer to https://github.com/Aseamann/ACEdecoy
import collections
import logging
import os
import pprint
import statistics
import sys
from math import sqrt
from pathlib import Path
from typing import Optional, List, Union, Tuple, final

import Bio
import numpy as np
from Bio import Align
from Bio.PDB import PDBParser
from scipy.spatial.transform import Rotation
from sklearn.decomposition import PCA
from sklearn.neighbors import KDTree

from ConfigManager import ConfigManager
from modules.wrappers.RosettaWrapper import REBprocessor
from util import file_utils

cm = ConfigManager.get_instance

TOOL_VER = "0.2.1_alpha"

# TODO: Make sure all methods accept Paths as well as str

# These types of PDB records will be carried / written over into the newly generated file when PDBs are modified
KEEP_LINES = ["HEADER", "REMARK", "SSBOND", "ATOM  ", "TER   ", "END   "]

# Magic numbers to extract data from PDB format, coordinate section, ATOM records.
# These specify the columns in which the relevant data is encoded in standard PDB files.
# See: https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
section_atom_number = [6, 11]
section_atom_name = [13, 16]
section_altloc = 16
section_residue = [17, 20]
section_chain_id = 21
section_residue_number = [22, 26]
section_x = [30, 38]
section_y = [38, 46]
section_z = [46, 54]
section_occupancy = [54, 60]
section_temp_factor = [60, 66]
section_element = [76, 78]
# Magic numbers to extract data from PDB format, connectivity section, SSBOND records.
# These specify the columns in which the relevant data is encoded in standard PDB files.
# See: https://www.wwpdb.org/documentation/file-format-content/format33/sect6.html#SSBOND
section_ssbond_residue_id1 = [17, 21]
section_ssbond_chain_id1 = 15
section_ssbond_residue_id2 = [31, 35]
section_ssbond_chain_id2 = 29


def get_chains(pdb: Union[str, Path]) -> Optional[list]:
    """
    Returns a list of all chains in PDB file.

    :param pdb: Path to PDB file to read
    :return: List of chains contained in PDB file, or None, if no chains can be found in PDB and VIPER is running in
        permissive mode
    """
    logging.debug(f"Getting chains from {pdb}...")
    chains = []
    with open(pdb, "r") as file:
        for line in file:
            if line[0:6] == "ATOM  ":
                if line[section_chain_id] not in chains:
                    chains.append(line[section_chain_id])
    logging.debug(f"Got chains '{chains}'")
    if len(chains) == 0:
        logging.error(f"Couldn't find any chain in '{pdb}'!")
        if cm().get("permissive"):
            return None
        else:
            sys.exit(1)
    return chains


def get_amino_acids_on_chain(pdb: Union[str, Path], chain: str, with_id: bool = False) -> Union[
    str, List[Tuple[str, int]]]:
    """
    Returns the string of amino acids in a specific chain as a string in single letter notation.

    :param pdb: Path to PDB file to read
    :param chain: Which chain to read
    :param with_id: Whether to also return the residue id for each amino acid.
    :return: String of amino acids in single letter formatting. The string may be empty if chain can't be found in PDB
        and VIPER is running in permissive mode. If with_id = True, instead returns a list of amino acid-residue id tuples
    """
    logging.debug(f"Trying to get amino acid sequence of chain '{chain}' in '{pdb}'...")
    output = ''
    if with_id:
        output = []
    flag = True
    with open(pdb, 'r') as file:
        for line in file:
            if line[0:6] == 'ATOM  ':
                if line[section_chain_id] == chain:
                    if flag:
                        count = int(line[section_residue_number[0]:section_residue_number[1]])
                        flag = False
                    if count == int(line[section_residue_number[0]:section_residue_number[1]]):
                        if line[section_altloc] != 'B':
                            if with_id:
                                output.append((three_to_one(line[section_residue[0]:section_residue[1]]), count))
                            else:
                                output += three_to_one(line[section_residue[0]:section_residue[1]])
                            count += 1
                    elif count < int(line[section_residue_number[0]:section_residue_number[1]]):
                        count = int(line[section_residue_number[0]:section_residue_number[1]])
    logging.debug(f"Got following sequence: '{output}'")
    if len(output) == 0:
        logging.error(f"Couldn't find any amino acid in chain '{chain}' in '{pdb}'!")
        if not cm().get("permissive"):
            raise ValueError(f"Couldn't find any amino acid in chain '{chain}' in '{pdb}'!")
    return output


def three_to_one(three: str, default: str = "") -> str:
    """
    Converts three letter amino acid abbreviation to single letter abbreviation.

    :param three: Three letter amino acid abbreviation
    :param default: Which optional string to return if the three letter abbreviation is not recognized
    :return: Single letter amino acid abbreviation, or '' or optional default if input abbreviation is not recognized
        and VIPER is running in permissive mode.
    """
    translate = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'ASX': 'B', 'CYS': 'C', 'GLU': 'E',
        'GLN': 'Q', 'GLX': 'Z', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y',
        'VAL': 'V'
    }
    if three.upper() in translate:
        return translate[three.upper()]
    elif three.upper() in translate.values():
        return three.upper()
    else:
        logging.log(logging.WARN, f"Tried to convert '{three}' to single letter amino acid abbreviation but failed!")
        if default:
            return default
        else:
            return ""


def one_to_three(one: str, default: str = "   ") -> str:
    """
    Converts three letter amino acid to single letter abbreviation.

    :param one: One letter amino acid abbreviation
    :param default: Which optional string to return if the one letter abbreviation is not recognized
    :return: Three letter amino acid abbreviation, or '   ' or optional default if the one letter amino acid
        abbreviation is not recognized, and VIPER is run in permissive mode.
    """
    translate = {
        'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'B': 'ASX', 'C': 'CYS', 'E': 'GLU',
        'Q': 'GLN', 'Z': 'GLX', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
        'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR',
        'V': 'VAL'
    }
    if one.upper() in translate:
        return translate[one.upper()]
    elif one.upper() in translate.values():
        return one.upper()
    else:
        logging.log(logging.WARN, f"Tried to convert '{one}' to three letter amino acid abbreviation but failed!")
        if default:
            return default
        else:
            return "   "


def first_atom_on_chain(pdb: str, chain: str) -> Optional[dict]:
    """
    Returns an atom dictionary for the first atom of a chain.

    :param pdb: Path to PDB to be read
    :param chain: Which chain to read
    :return: A dict of properties  of the first atom in the specified chain, or None, if no first atom could be found
        and VIPER is running in permissive mode
    """
    logging.debug(f"Trying to get first atom on chain '{chain}' in '{pdb}'...")
    with open(pdb, 'r') as file:
        for line in file:
            if line[0:6] == 'ATOM  ':
                if line[section_chain_id] == chain.upper():
                    atom = {'atom_num': int(line[section_atom_number[0]:section_atom_number[1]]),
                            'atom_id': line[section_atom_name[0]:section_atom_number[1]].strip(),
                            'atom_comp_id': line[section_residue[0]:section_residue[1]],
                            'chain_id': line[section_chain_id],
                            'comp_num': int(line[section_residue_number[0]:section_residue_number[1]]),
                            'X': float(line[section_x[0]:section_x[1]]),
                            'Y': float(line[section_y[0]:section_y[1]]),
                            'Z': float(line[section_z[0]:section_z[1]]),
                            'occupancy': float(line[section_occupancy[0]:section_occupancy[1]])}
                    if len(line) >= 76:
                        atom['B_iso_or_equiv'] = float(line[section_temp_factor[0]:section_temp_factor[1]])
                        atom['atom_type'] = line[section_element[0]:section_element[1]]
                    logging.log(logging.DEBUG, f"Found following atom: {atom}")
                    return atom
        logging.error(f"Couldn't find first atom on chain '{chain}' in '{pdb}'!")
        if cm().get("permissive"):
            return None
        else:
            sys.exit(1)


def get_atom(pdb: str, atom_num: int) -> Optional[dict]:
    """
    Returns a dictionary with elements related to specific atom number entered.

    :param pdb: Path to PDB file to be read
    :param atom_num: Integer value for atom in file (serial number)
    :return: Dictionary containing information for specified atom from PDB file, or None, if the atom could not be found
        and VIPER is running in permissive mode
    """
    logging.log(logging.DEBUG, f"Trying to get atom number {atom_num} from '{pdb}'...")
    with open(pdb, 'r') as file:
        for line in file:
            if line[0:6] == 'ATOM  ':
                if int(line[section_atom_number[0]:section_atom_number[1]]) == atom_num and len(line) >= 76:
                    atom = {'atom_num': int(line[section_atom_number[0]:section_atom_number[1]]),
                            'atom_id': line[section_atom_name[0]:section_atom_name[1]].strip(),
                            'atom_comp_id': line[section_residue[0]:section_residue[1]],
                            'chain_id': line[section_chain_id],
                            'comp_num': int(line[section_residue_number[0]:section_residue_number[1]]),
                            'X': float(line[section_x[0]:section_x[1]]),
                            'Y': float(line[section_y[0]:section_y[1]]),
                            'Z': float(line[section_z[0]:section_z[1]]),
                            'occupancy': float(line[section_occupancy[0]:section_occupancy[1]])}
                    if len(line) >= 76:
                        atom['B_iso_or_equiv'] = float(line[section_temp_factor[0]:section_temp_factor[1]])
                        atom['atom_type'] = line[section_element[0]:section_element[1]]
                    logging.log(logging.DEBUG, f"Found following atom: {atom}")
                    return atom
        logging.error(f"Couldn't locate atom number {atom_num} from '{pdb}'!")
        if cm().get("permissive"):
            return None
        else:
            sys.exit(1)


def get_atoms_on_chain(pdb: Union[str, Path], chain: str) -> list:
    """
    Collect the atoms from an inputted chain. Provides all values that are in the specified PDB file.

    :param pdb: Path to PDB to be read
    :param chain: Which chain to use
    :return: List of atoms based on the chain id that was submitted. This list may be empty if no atoms could be found
        and VIPER is running permissive mode
    """
    logging.debug(f"Trying to get atoms on chain '{chain}' in '{pdb}'...")
    atoms = []
    with open(os.path.normpath(pdb), 'r') as file:
        for line in file:
            if line[0:6] == 'ATOM  ':
                if line[section_chain_id] == chain and len(line) >= 76:
                    atom = {'atom_num': int(line[section_atom_number[0]:section_atom_number[1]]),
                            'atom_id': line[section_atom_name[0]:section_atom_name[1]].strip(),
                            'atom_comp_id': line[section_residue[0]:section_residue[1]],
                            'chain_id': line[section_chain_id],
                            'comp_num': int(line[section_residue_number[0]:section_residue_number[1]]),
                            'X': float(line[section_x[0]:section_x[1]]),
                            'Y': float(line[section_y[0]:section_y[1]]),
                            'Z': float(line[section_z[0]:section_z[1]]),
                            'occupancy': float(line[section_occupancy[0]:section_occupancy[1]])}
                    if len(line) >= 76:
                        atom['B_iso_or_equiv'] = float(line[section_temp_factor[0]:section_temp_factor[1]])
                        atom['atom_type'] = line[section_element[0]:section_element[1]]
                    atoms.append(atom)
    if len(atoms) == 0:
        logging.error(f"Tried to get atoms on chain '{chain}' from '{pdb}' but failed!")
        if not cm().get("permissive"):
            sys.exit(1)
    else:
        logging.debug(f"Atoms on chain '{chain}' are {atoms}.")
    return atoms


def euclidean_of_atoms(pdb: str, atom_num_1: int, atom_num_2: int) -> Optional[float]:
    """
    Returns the Euclidean distance between two atoms identified through their (serial) number.

    :param pdb: Path to PDB to be read
    :param atom_num_1: Atom number of first atom
    :param atom_num_2: Atom number of second atom
    :return: Distance between two atoms in angstroms, or None, if the distance could not be computed and VIPER is
        running in permissive mode
    """
    logging.debug(f"Trying to find euclidean distance between atom {atom_num_1} and {atom_num_2} in '{pdb}'...")
    atom_1 = get_atom(pdb, atom_num_1)
    atom_2 = get_atom(pdb, atom_num_2)
    if atom_1 is None or atom_2 is None:
        logging.error(
            f"Tried to compute euclidean distance between atom {atom_num_1} and {atom_num_2} in '{pdb}' but failed!"
            f"One of the two atoms could not be found!")
        if cm().get("permissive"):
            return None
        else:
            raise ValueError(
                f"Tried to compute euclidean distance between atom {atom_num_1} and {atom_num_2} in '{pdb}' but failed!"
                f"One of the two atoms could not be found!")
    else:
        euclidean_distance = sqrt(
            (atom_2['X'] - atom_1['X']) ** 2 + (atom_2['Y'] - atom_1['Y']) ** 2 + (atom_2['Z'] - atom_1['Z']) ** 2)
        logging.log(logging.DEBUG, f"Euclidean distance (Angstrom): {euclidean_distance}")
        return euclidean_distance


def rebuild_atom_line(atoms: list) -> str:
    """
    Rebuild the lines of atoms in a PDB based on their atom dictionaries

    :param atoms: List atom-type dicts (see get_atom())
    :return: Lines in a PDB corresponding to passed atom
    """
    if len(atoms) == 0:
        if cm().get("permissive"):
            return ""
        else:
            raise ValueError("Can't rebuild empty list!")
    last_chain = atoms[0]["chain_id"]
    output = ""
    for atom in atoms:
        line = 'ATOM  ' + str(atom['atom_num']).rjust(5) + "  " + atom['atom_id'].ljust(3) + " "
        line += atom['atom_comp_id'] + " " + atom['chain_id'] + str(atom['comp_num']).rjust(4) + "    "
        line += str(format(atom['X'], ".3f")).rjust(8) + str(format(atom['Y'], ".3f")).rjust(8)
        line += str(format(atom['Z'], ".3f")).rjust(8)
        line += str(format(atom['occupancy'], ".2f")).rjust(6)
        line += str(format(atom['B_iso_or_equiv'], ".2f")).rjust(6) + "           "
        line += atom['atom_type'] + "\n"
        output += line
        if last_chain != atom["chain_id"]:
            last_chain = atom["chain_id"]
            output += "TER\n"
    logging.debug(f"Rebuilt atom(s): {output}")
    return output


def renumber_ascending(pdb: str, out_path: str = None) -> Path:
    """
    Write out a PDB file in the same directory as the input PDB with atom number and sequence number in ascending order.
    There is an inherent limit on how large atom (serial) ids and residue ids may be (99,999 and 9,999).
    If the proteins are too large, the count may exceed this number. Should VIPER be running in permissive mode, the ids
    will be taken modulo 100,000 and 10,000, respectively.
    This may cause errors later on, since ids will be duplicated!
    Note: Only HEADER, SSBOND, and ATOM records will be kept, and no secondary positions will be copied over

    :param pdb: Path to PDB to be read
    :param out_path: Optional path for PDB file being created, otherwise "_renumbered" appended to original name and the
        PDB will be written to the same directory
    :return: A Path object to the written out file
    """
    logging.info(f"Renumbering {pdb} by counting up...")

    out_name = pdb[:-4] + "_renumbered.pdb"
    if out_path:
        out_name = out_path

    header = (
        f"EXPDTA     MODEL                  RENUMBERED\nREMARK    3                                 \nREMARK    3 "
        f"PROGRAM    : VIPER, PDBtool {TOOL_VER} \n")
    ssbond = ""
    ssbond_res_ids = {}  # Will become mapping for updated ids
    coordinates = ""

    running_atom_count = 1  # Keeps track of atom number
    previous_residue_number = -10000
    running_residue_count = 0  # Keeps track of residue number
    chains = []  # Chains already seen
    with open(pdb, "r") as i:
        for line in i:
            if line[0:6] == 'HEADER':
                header = line + header
            if line[0:6] == "SSBOND":
                ssbond += line
                ssbond_res_ids[(int(line[section_ssbond_residue_id1[0]:section_ssbond_residue_id1[1]]),
                                line[section_ssbond_chain_id1])] = -1
                ssbond_res_ids[(int(line[section_ssbond_residue_id2[0]:section_ssbond_residue_id2[1]]),
                                line[section_ssbond_chain_id2])] = -1
            if line[0:6] == 'ATOM  ':  # Only copy over atoms
                if line[16] != 'B' and line[26] == ' ':  # Don't allow secondary atoms
                    current_residue_number = int(line[section_residue_number[0]:section_residue_number[1]])
                    if len(chains) == 0:
                        chains.append(line[section_chain_id])
                    elif line[section_chain_id] != chains[-1]:  # Encountered new chain
                        chains.append(line[section_chain_id])
                        coordinates += "TER\n"
                    if current_residue_number != previous_residue_number:  # Encountered next residue
                        previous_residue_number = int(line[section_residue_number[0]:section_residue_number[1]])
                        running_residue_count += 1
                    # Update SSBOND mapping with new res_id
                    for k in ssbond_res_ids.keys():
                        if current_residue_number == k[0] and line[section_chain_id] == k[1]:
                            ssbond_res_ids[(current_residue_number, line[section_chain_id])] = running_residue_count
                    # Update atom count and residue count
                    # Because there are only 5 characters for the atom serial number and 4 for the residue, we need
                    # to restrict the number space to < 100,000 and 10,000, respectively
                    if running_atom_count > 99999 or running_residue_count > 9999:
                        if cm().get("permissive"):
                            logging.log(logging.WARN, f"While renumbering {pdb} atom or residue count exceeded legal "
                                                      f"limits. Counts will be taken modulo 100,000 and 10,000, "
                                                      f"respectively. Proceeding may cause unexpected program "
                                                      f"behaviour.")
                            line = line[:6] + str(running_atom_count % 100000).rjust(5) + line[11:22] + str(
                                running_residue_count % 10000).rjust(4) + line[26:]
                        else:
                            logging.log(logging.ERROR, f"While renumbering {pdb} atom or residue count exceeded legal "
                                                       f"limits (100000 or 10000, respectively. Aborting! You might "
                                                       f"want to try to reduce the proteins to a cut out version "
                                                       f"focusing on the binding site.")
                            sys.exit(1)
                    else:
                        line = line[:6] + str(running_atom_count).rjust(5) + line[11:22] \
                               + str(running_residue_count).rjust(4) + line[26:]

                    running_atom_count += 1
                    coordinates += line
        coordinates += "TER\nEND\n"
    for old_entry, new_entry in ssbond_res_ids.items():
        old = old_entry[1] + " " + str(old_entry[0]).rjust(4)
        new = old_entry[1] + " " + str(new_entry).rjust(4)
        ssbond = ssbond.replace(old, new, 1)
    with open(out_name, "w+") as o:
        o.write(header)
        o.write(ssbond)
        o.write(coordinates)
    return Path(out_name)


def remove_chain(pdb: str, chain_id: List[str], rename: str = None) -> None:
    """
    Removes chains with a certain ID and save the output in a PDB.

    :param pdb: Path to PDB to be read
    :param chain_id: A list of ids to be deleted (i.e. ['A', 'E', 'P'])
    :param rename: Optional name for PDB file being created, otherwise use "_removed_chains<ids>" appended to original name
    """
    ids = ''.join(chain_id).upper()
    logging.info(f"Trying to remove chains {ids} from {pdb}...")
    out_name = pdb[:-4] + f"_removed_chains{ids}.pdb"
    if rename:
        out_name = rename
    with open(pdb, 'r') as i, open(out_name, 'w+') as o:
        for line in i:
            if line[0:6] in KEEP_LINES:
                if line[0:6] == "ATOM  " and line[section_chain_id].upper() not in ids:
                    o.write(line)
                    continue
                if (line[0:6] == "SSBOND" and line[section_ssbond_chain_id1] not in ids
                        and line[section_ssbond_chain_id2] not in ids):
                    o.write(line)
                    continue
                o.write(line)  # is a line where chain id doesn't matter (HEADER, REMARK, ...)
        logging.info(f"PDB with chains {ids} removed saved to {rename}!")


def superimpose_multiple(pdb: str, ref_pdb: str, target_order: str, ref_order: str,
                         path: Union[str, Path, List[Union[str, Path]]]) -> Tuple[Path, float]:
    """
    Superimpose two PDBs.
    Read in target and reference structure, superimpose target to reference, save superimposed target structure.
    Writes the result to 'path' appended to the results path, creating any subdirectories as needed.
    The path will be based off the results path, so do not use absolute paths!
    Instead, use a list with directory names and end it with the filename.

    :param pdb: Location of target PDB
    :param ref_pdb: Location of reference PDB
    :param target_order: Order of chains to compare in superimposing structures for target
    :param ref_order: Order of chains to compare in superimposing structures for reference
    :param path: Where to write the file to
    :return RMSD value for all-atom RMSD
    """
    raise DeprecationWarning
    # This was for testing and will be removed in the future
    aligner = Bio.Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.gap_score = -100.00
    aligner.match_score = 2.0
    aligner.mismatch_score = 0.0

    logging.warning(f"Please use superimpose_single() multiple times instead!")

    logging.info(f"Trying to superimpose {pdb} onto {ref_pdb} with chain order {target_order} (target) and "
                 f"{ref_order} (reference).")

    # Collect target pdb seq
    target_seq = {}  # Dictionary chain_id:seq
    for chain in get_chains(pdb):
        residues = get_amino_acids_on_chain(pdb, chain)
        if residues is None:
            logging.error(f"Couldn't find any residues for chain {chain} in '{pdb}'!")
            if not cm().get("permissive"):
                logging.error(f"Aborting now...")
                sys.exit(1)
        target_seq[chain] = residues

    # Collect reference pdb seq
    ref_seq = {}  # Dictionary chain_id:seq
    for chain in get_chains(ref_pdb):
        residues = get_amino_acids_on_chain(ref_pdb, chain)
        if residues is None:
            logging.error(f"Couldn't find any residues for chain {chain} in '{pdb}'!")
            if not cm().get("permissive"):
                logging.error(f"Aborting now...")
                sys.exit(1)
        ref_seq[chain] = residues

    # Starting positions of each amino acid from target and reference
    # 'target' -> 'chain' -> [start, end]
    start_pos = {"target": {}, "reference": {}}
    for chain in list(target_order):
        # print(chain)
        start_pos['target'][chain] = []
    for chain in list(ref_order):
        # print(chain)
        start_pos['reference'][chain] = []

    # Produce alignments
    target_pos = list(target_order)
    ref_pos = list(ref_order)
    for chain_pos in range(0, len(target_pos)):
        # Align with chain in target_order and ref_order
        target_chain = target_seq[target_pos[chain_pos]]
        reference_chain = ref_seq[ref_pos[chain_pos]]
        temp_align = aligner.align(target_chain, reference_chain)
        logging.debug(f"Target Chain: {target_pos[chain_pos]}")
        logging.debug(f"Reference Chain: {ref_pos[chain_pos]}")
        logging.debug(f"{temp_align[0]}")

        start_pos['target'][target_pos[chain_pos]] = [temp_align[0].path[0][0], temp_align[0].path[1][0]]
        start_pos['reference'][ref_pos[chain_pos]] = [temp_align[0].path[0][1], temp_align[0].path[1][1]]

    # Initialize parser
    if cm().get("verbose"):
        parser = Bio.PDB.PDBParser(QUIET=False)
    else:
        parser = Bio.PDB.PDBParser(QUIET=True)

    # Gather Structures
    ref_structure = parser.get_structure("reference", ref_pdb)
    target_structure = parser.get_structure("target", pdb)

    # Collect structures
    ref_model = ref_structure[0]
    target_model = target_structure[0]

    # Collecting alpha carbons
    ref_atoms = []
    target_atoms = []
    # Reference Model
    for ref_chain in ref_model:
        chain = ref_chain.__repr__().split("=")[1].split(">")[0]
        # print(ref_chain)
        # print(chain)
        if chain in ref_pos:
            first_in_chain = True
            skip_switch = False
            start = 0
            end = 0
            last_count = 0
            for ref_res in ref_chain:
                if ref_res.get_resname() != "HOH" and 'CA' in ref_res:
                    if first_in_chain:
                        offset = ref_res.get_id()[1]
                        start = start_pos['reference'][chain][0] + offset
                        end = start_pos['reference'][chain][1] + offset
                        first_in_chain = False
                    if ref_res.get_id()[1] in range(start, end):
                        if skip_switch:
                            if ref_res.get_id()[1] != last_count + 1:  # Catch when count skips a number
                                end += ref_res.get_id()[1] - last_count - 1
                        ref_atoms.append(ref_res['CA'])
                        last_count = ref_res.get_id()[1]
                        skip_switch = True
    # Target Model
    for target_chain in target_model:
        chain = target_chain.__repr__().split("=")[1].split(">")[0]
        if chain in target_pos:
            first_in_chain = True
            skip_switch = False
            start = 0
            end = 0
            last_count = 0
            for target_res in target_chain:
                if target_res.get_resname() != "HOH" and 'CA' in target_res:
                    if first_in_chain:  # If chain doesn't start count at position 0_old
                        offset = target_res.get_id()[1]
                        start = start_pos['target'][chain][0] + offset
                        end = start_pos['target'][chain][1] + offset
                        first_in_chain = False
                    if target_res.get_id()[1] in range(start, end):
                        if skip_switch:
                            if target_res.get_id()[1] != last_count + 1:  # Catch when count skips a number
                                end += target_res.get_id()[1] - last_count - 1
                        target_atoms.append(target_res['CA'])
                        last_count = target_res.get_id()[1]
                        skip_switch = True
    # print(len(ref_atoms))
    # print(len(target_atoms))
    # Truncate long list
    if len(ref_atoms) > len(target_atoms):
        ref_atoms = ref_atoms[:len(target_atoms)]
    elif len(target_atoms) > len(ref_atoms):
        target_atoms = target_atoms[:len(ref_atoms)]
    # Superimposing
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, target_atoms)
    super_imposer.apply(target_model.get_atoms())

    # Save structure
    io = Bio.PDB.PDBIO()
    logging.info(
        f"Stats for superimposed structures - RMSD: {super_imposer.rms}, Atoms Pulled: {str(len(target_atoms))}")
    io.set_structure(target_structure)

    p = file_utils.make_file(path, "", dry_run=True)
    io.save(p)
    return p, super_imposer.rms


def superimpose_single(pdb: Union[str, Path], ref_pdb: Union[str, Path], query_chain: str, ref_chain: str,
                       out_path: Union[str, Path, List[Union[str, Path]]],
                       aligner: Bio.Align.PairwiseAligner = None) -> Tuple[Path, float]:
    """
    Superimposes a single chain from pdb on to a single chain from ref_pdb. Performs an alignment first, such that the
    largest possible subsequence without creating a gap while accepting mismatches gets aligned, unless a custom
    aligner is supplied. When using the standard, this only uses the alpha carbons to superimpose the backbone.
    Writes only the superimposed pdb chain to out_path, a file/list of directories + file in the overall results_path.

    :param pdb: The pdb from which to take the chain that should be superimposed
    :param ref_pdb: The reference pdb, on which the chain should be superimposed
    :param query_chain: Which chain from pdb to superimpose
    :param ref_chain: Which chain from ref_pdb to superimpose on to
    :param out_path: Where to write out the superimposed chain PDB file
    :param aligner: (Optional) A custom aligner for determining which atoms should be used for calculating the
        translation and rotation to superimpose on to the reference
    :return: A tuple of the path to where the PDB was written and the superimposer RMSD
    """
    logging.info(f"Trying to superimpose {pdb}, chain {query_chain} onto {ref_pdb}, chain {ref_chain}...")
    try:
        query_pdb = Path(pdb).resolve(strict=True)
        ref_pdb = Path(ref_pdb).resolve(strict=True)
    except (FileNotFoundError, RuntimeError) as e:
        logging.error(f"Couldn't resolve a path you provided. Stacktrace: {e}")
        raise e

    pdb_parser = Bio.PDB.PDBParser()
    query_chain_struc = pdb_parser.get_structure(query_pdb.name, query_pdb)[0][query_chain]
    ref_chain_struc = pdb_parser.get_structure(ref_pdb.name, ref_pdb)[0][ref_chain]
    query_seq = [(residue.id, three_to_one(residue.resname, "X")) for residue in query_chain_struc.get_residues()]
    ref_seq = [(residue.id, three_to_one(residue.resname, "X")) for residue in ref_chain_struc.get_residues()]
    query_aaseq = "".join([residue[1] for residue in query_seq])
    ref_aaseq = "".join([residue[1] for residue in ref_seq])

    if aligner is None:
        aligner = Bio.Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.gap_score = -2.0
        aligner.match_score = 2.0
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -.5
    alignment = aligner.align(ref_aaseq, query_aaseq)  # where in ref is target?
    # Determine best alignment
    curr_best = None
    for a in alignment:
        if curr_best is None:
            curr_best = a
            continue
        if a.score > curr_best.score:
            curr_best = a
    logging.debug(f"Best alignment: {str(curr_best)}")

    query_aligned = curr_best.query  # We're querying the ref chain as to where the target chain is, so that we know at
    ref_aligned = curr_best.target  # which residue in ref the target chain appears, therefore _target_ is the reference
    # Which residue id in the reference corresponds to which residue id in the queried pdb?
    residue_id_map = {}
    # .aligned has format:
    # [[[ target_start, target_end ]] [[ query_start, query_end ]]]
    query_start, query_end = curr_best.aligned[1][0][0], curr_best.aligned[1][0][1]
    ref_start, ref_end = curr_best.aligned[0][0][0], curr_best.aligned[0][0][1]

    query_offset = query_start
    ref_offset = ref_start
    for (ref_aa, query_aa) in zip((ref_aligned[ref_start:ref_end]), query_aligned[query_start:query_end]):
        # Since the penalty for gaps is so high, there probably shouldn't be any in the output. This is just for
        # safety
        if ref_aa == "-":
            if query_aa != "-":  # Gap in ref
                query_offset += 1
                continue
        if query_aa == "-":
            if ref_aa != "-":  # Gap in query
                ref_offset += 1
                continue
        # We're allowing mismatches within the sequence to enable us to modify the peptide later on and still align it
        # Therefore we might map different amino acids on to each other!
        residue_id_map[ref_seq[ref_offset][0]] = query_seq[query_offset][0]
        query_offset += 1
        ref_offset += 1
    logging.debug(f"The residue id mapping (ref -> query) is: {pprint.pformat(residue_id_map)}")

    ref_atom_list = []
    query_atom_list = []
    for residue in residue_id_map:
        # Only add alpha carbons to align backbone, maybe relaxation is necessary afterwards?
        ref_atom_list.append(ref_chain_struc[residue]["CA"])
        query_atom_list.append(query_chain_struc[residue_id_map[residue]]["CA"])

    superimposer = Bio.PDB.Superimposer()
    superimposer.set_atoms(ref_atom_list, query_atom_list)
    superimposer.apply(query_chain_struc.get_atoms())
    logging.info(f"Successfully superimposed structures with RMSD: {superimposer.rms:+.4f}")

    p = file_utils.make_file(out_path, "", dry_run=True)
    logging.info(f"Writing out superimposed pdb to '{p}'")
    pdb_io = Bio.PDB.PDBIO()
    pdb_io.set_structure(query_chain_struc)
    pdb_io.save(os.path.normpath(p))

    return p, superimposer.rms


def superimpose_reflist(nodes: List[REBprocessor.Node], query_pdb: Union[str, Path], ref_pdb: Union[str, Path],
                        query_chain: str, ref_chain: str,
                        out_path: Union[str, Path, List[Union[str, Path]]],
                        use_longest_subseq: bool = True) -> Tuple[Path, float]:
    """
    Tries to superimpose two PDBs using a reference list of nodes. Nodes hold an orig_res_id field, which indicates
    what nodes they originally were when a PDB was read in. This method assumes that the nodes in 'nodes' are in the
    same order as the residues in 'query_pdb'. It then tries to superimpose the alpha carbons of all residues from
    'query_pdb' onto 'ref_pdb' by using the reference 'nodes' to link the residues in 'query_pdb' to 'ref_pdb'.
    Specifically, it calculates the singular value decomposition to transform the alpha carbons from 'query_pdb' onto
    the alpha carbons in 'ref_pdb' that correspond to each other. If 'use_longest_subseqs' is True (default), it only
    uses the alpha carbons of the longest continuous subsequence in 'nodes' to calculate the singular value
    decomposition, otherwise it uses all alpha carbons it can match to 'ref_pdb'.

    :param nodes: A reference list of REBprocessor.Node objects to map query onto ref
    :param query_pdb: The PDB to superimpose onto the reference
    :param ref_pdb: The reference to superimpose onto
    :param query_chain: Which chain in query_pdb to use
    :param ref_chain: Which chain in ref_pdb to use
    :param out_path: Where to save the resultant superimposed query_pdb to (this gets appended to the path to the
        general VIPER output folder)
    :param use_longest_subseq: Whether to only use the longest continuous subsequence in nodes to calculate the singular
        value decomposition (default: True)
    :return: A tuple of the path to the written out pdb and the superimposer RMSD
    """
    raise DeprecationWarning
    # This was for testing and will be removed in the future
    logging.info(f"Trying to superimpose {query_pdb}, chain {query_chain} onto {ref_pdb}, chain {ref_chain}...")
    try:
        query_pdb = Path(query_pdb).resolve(strict=True)
        ref_pdb = Path(ref_pdb).resolve(strict=True)
    except (FileNotFoundError, RuntimeError) as e:
        logging.error(f"Couldn't resolve a path you provided. Stacktrace: {e}")
        raise e

    query_atoms = get_atoms_on_chain(query_pdb, query_chain)
    ref_atoms = get_atoms_on_chain(ref_pdb, ref_chain)

    atom_positions = {}
    use_nodes = [(n, node) for (n, node) in enumerate(nodes, start=1)]

    # Determine the longest continuous subsequence of nodes
    if use_longest_subseq:
        subseqs = [[]]
        offset = 0
        first = True
        # Forward scan
        for (n, node) in use_nodes:
            if first:
                subseqs[offset].append((n, node))
                first = False
                continue
            if node.orig_res_id is None:
                continue
            orig_id = node.orig_res_id[0]
            if orig_id - subseqs[offset][-1][1].orig_res_id[0] > 1:
                # Encountered gap in sequence, next node's id is more than 1 away. Start new subsequence
                subseqs.append([])
                offset += 1
            subseqs[offset].append((n, node))
        use_nodes = max(subseqs, key=lambda l: len(l))  # restrict to longest subsequence

    # Extract coordinates of alpha carbons
    for n, node in use_nodes:
        if node.orig_res_id is not None:
            nid = node.orig_res_id
            if nid[1] != query_chain:
                raise ValueError(
                    f"The chain id of passed node {node} is not the same as the passed query_chain '{query_chain}'")
            if not atom_positions.get(nid[0], False):
                atom_positions[nid[0]] = {}
            for atom in query_atoms:
                if atom["comp_num"] == n and atom["atom_id"] == "CA":
                    atom_positions[nid[0]]["query"] = [atom["X"], atom["Y"], atom["Z"]]
            for atom in ref_atoms:
                if atom["comp_num"] == nid[0] and atom["atom_id"] == "CA":
                    atom_positions[nid[0]]["ref"] = [atom["X"], atom["Y"], atom["Z"]]

    # Match query alpha carbons to reference alpha carbons
    query_pos_list = []
    ref_pos_list = []
    for id, positions in atom_positions.items():
        q = positions.get("query")
        r = positions.get("ref")
        if q is not None and r is not None:
            query_pos_list.append(q)
            ref_pos_list.append(r)

    query = np.array(query_pos_list, "f")
    ref = np.array(ref_pos_list, "f")

    # Use BioPython singular value decomposition wrapper
    svds = Bio.SVDSuperimposer.SVDSuperimposer()
    svds.set(ref, query)
    svds.run()
    rms = svds.get_rms()
    rotation, translation = svds.get_rotran()
    for atom in query_atoms:
        coords = [[atom["X"], atom["Y"], atom["Z"]]]
        coords = np.dot(np.array(coords), rotation) + translation
        atom["X"] = float(coords[0][0])
        atom["Y"] = float(coords[0][1])
        atom["Z"] = float(coords[0][2])

    p = file_utils.make_file(path=out_path, content=rebuild_atom_line(query_atoms))
    return p, rms


def kabsch(nodes: List[REBprocessor.Node], query_pdb: Union[str, Path], ref_pdb: Union[str, Path],
           query_chain: str, ref_chain: str,
           out_path: Union[str, Path, List[Union[str, Path]]],
           use_longest_subseq: bool = True) -> Tuple[Path, float]:
    """
    Tries to superimpose two PDBs using a reference list of nodes. Nodes hold an orig_res_id field, which indicates
    what nodes they originally were when a PDB was read in. This method assumes that the nodes in 'nodes' are in the
    same order as the residues in 'query_pdb'. It then tries to superimpose the alpha carbons of all residues from
    'query_pdb' onto 'ref_pdb' by using the reference 'nodes' to link the residues in 'query_pdb' to 'ref_pdb'.
    Specifically, it uses the Kabsch algorithm to transform the alpha carbons from 'query_pdb' onto the alpha carbons
    in 'ref_pdb' that correspond to each other. If 'use_longest_subseqs' is True (default), it only uses the alpha
    carbons of the longest continuous subsequence in 'nodes' to calculate the singular value decomposition, otherwise
    it uses all alpha carbons it can match to 'ref_pdb'.

    :param nodes: A reference list of REBprocessor.Node objects to map query onto ref
    :param query_pdb: The PDB to superimpose onto the reference
    :param ref_pdb: The reference to superimpose onto
    :param query_chain: Which chain in query_pdb to use
    :param ref_chain: Which chain in ref_pdb to use
    :param out_path: Where to save the resultant superimposed query_pdb to (this gets appended to the path to the
        general VIPER output folder)
    :param use_longest_subseq: Whether to only use the longest continuous subsequence in nodes to calculate the singular
        value decomposition (default: True)
    :return: A tuple of the path to the written out pdb and the superimposer RMSD
    """
    raise DeprecationWarning
    # This was for testing and will be removed in the future
    logging.info(f"Trying to superimpose {query_pdb}, chain {query_chain} onto {ref_pdb}, chain {ref_chain}...")
    try:
        query_pdb = Path(query_pdb).resolve(strict=True)
        ref_pdb = Path(ref_pdb).resolve(strict=True)
    except (FileNotFoundError, RuntimeError) as e:
        logging.error(f"Couldn't resolve a path you provided. Stacktrace: {e}")
        raise e

    query_atoms = get_atoms_on_chain(query_pdb, query_chain)
    ref_atoms = get_atoms_on_chain(ref_pdb, ref_chain)

    atom_positions = {}

    use_nodes = [(n, node) for (n, node) in enumerate(nodes, start=1)]

    # Determine the longest continuous subsequence of nodes
    if use_longest_subseq:
        subseqs = [[]]
        offset = 0
        first = True
        # Forward scan
        for (n, node) in use_nodes:
            if first:
                subseqs[offset].append((n, node))
                first = False
                continue
            if node.orig_res_id is None:
                continue
            orig_id = node.orig_res_id[0]
            if orig_id - subseqs[offset][-1][1].orig_res_id[0] > 1:
                # Encountered gap in sequence, next node's id is more than 1 away. Start new subsequence
                subseqs.append([])
                offset += 1
            subseqs[offset].append((n, node))
        use_nodes = max(subseqs, key=lambda l: len(l))  # restrict to longest subsequence

    # Extract coordinates of alpha carbons
    for n, node in use_nodes:
        if node.orig_res_id is not None:
            nid = node.orig_res_id
            if nid[1] != query_chain:
                raise ValueError(
                    f"The chain id of passed node {node} is not the same as the passed query_chain '{query_chain}'")
            if not atom_positions.get(nid[0], False):
                atom_positions[nid[0]] = {}
            for atom in query_atoms:
                if atom["comp_num"] == n and atom["atom_id"] == "CA":
                    atom_positions[nid[0]]["query"] = [atom["X"], atom["Y"], atom["Z"]]
            for atom in ref_atoms:
                if atom["comp_num"] == nid[0] and atom["atom_id"] == "CA":
                    atom_positions[nid[0]]["ref"] = [atom["X"], atom["Y"], atom["Z"]]

    # Match query alpha carbons to reference alpha carbons
    query_pos_list = []
    ref_pos_list = []
    for id, positions in atom_positions.items():
        q = positions.get("query")
        r = positions.get("ref")
        if q is not None and r is not None:
            query_pos_list.append(q)
            ref_pos_list.append(r)

    query = np.array(query_pos_list, "f")
    ref = np.array(ref_pos_list, "f")
    q_mean = np.mean(query, axis=0)
    r_mean = np.mean(ref, axis=0)

    # Center the data points
    q_centered = query - q_mean
    r_centered = ref - r_mean

    # Calculate the covariance matrix
    H = np.dot(q_centered.T, r_centered)

    # Singular Value Decomposition
    U, _, Vt = np.linalg.svd(H)

    # Ensure proper rotation matrix (no reflection) with correct determinant
    D = np.eye(3)
    D[2, 2] = np.linalg.det(np.dot(U, Vt))

    # Calculate the rotation matrix
    R = np.dot(U, np.dot(D, Vt))

    # Only applying the rotation (centered on origin)
    for atom in query_atoms:
        coords = [[atom["X"], atom["Y"], atom["Z"]]]
        coords = np.dot(np.array(coords), R)
        atom["X"] = float(coords[0][0])
        atom["Y"] = float(coords[0][1])
        atom["Z"] = float(coords[0][2])

    p = file_utils.make_file(path=out_path, content=rebuild_atom_line(query_atoms))

    # Add the mean of the reference positions to translate to original position
    for atom in query_atoms:
        coords = [[atom["X"], atom["Y"], atom["Z"]]]
        coords = np.array(coords) + r_mean
        atom["X"] = float(coords[0][0])
        atom["Y"] = float(coords[0][1])
        atom["Z"] = float(coords[0][2])

    out_path[-1] = out_path[-1][:-4] + "_addmean.pdb"
    p = file_utils.make_file(path=out_path, content=rebuild_atom_line(query_atoms))

    # Subtract the rotated query mean to offset accordingly
    for atom in query_atoms:
        coords = [[atom["X"], atom["Y"], atom["Z"]]]
        coords = np.array(coords) - np.dot(q_mean, R)
        atom["X"] = float(coords[0][0])
        atom["Y"] = float(coords[0][1])
        atom["Z"] = float(coords[0][2])

    out_path[-1] = out_path[-1][:-4] + "_subquery.pdb"
    p = file_utils.make_file(path=out_path, content=rebuild_atom_line(query_atoms))

    # TODO: Add RMSE calculation
    return Path("."), 0.0


def rmsd(pdb: str, ref_pdb: str, target_order: str, ref_order: str, ca: bool = False) -> float:
    """
    Calculate RMSD values between two PDBs - based on aligned amino acids. Optional Carbon Alpha RMSD as well.

    :param pdb: Path to target PDB
    :param ref_pdb: Reference PDB - can be same PDB as target
    :param target_order: Target chains to consider for RMSD calculation - must be in same order as reference
    :param ref_order: Reference chains to consider for RMSD calculation - must be in same order as target
    :param ca: True - only Alpha Carbon RMSD ; False - all-atom RMSD
    :return: RMSD value
    """

    logging.info(f"Calculating RMSD for '{pdb}' and '{ref_pdb}' with chains '{target_order}' and '{ref_order}'"
                 f", respectively...")

    target_atoms = {}  # chain: atoms
    ref_atoms = {}  # chain: atoms
    target_aa = {}  # chain: amino acids
    ref_aa = {}  # chain: amino acids
    rmsd_array = {"target": [], "ref": []}  # List of cords in order
    # Collect target atoms
    for chain in target_order:
        # Collect atoms - contains xyz and aa
        target_atoms[chain] = get_atoms_on_chain(pdb, chain)
        # Collect just residue list
        target_aa[chain] = get_amino_acids_on_chain(pdb, chain)
    # Collect reference atoms
    for chain in ref_order:
        # Collect atoms - contains xyz and aa
        ref_atoms[chain] = get_atoms_on_chain(ref_pdb, chain)
        # Collect just residue list
        ref_aa[chain] = get_amino_acids_on_chain(ref_pdb, chain)
    # Return to previous PDB
    # Run alignment
    aligner = Bio.Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.gap_score = -100.00
    aligner.match_score = 2.0
    aligner.mismatch_score = 0.0
    # Loop through each paired chain
    for position in range(len(target_order)):
        # run alignment
        temp_align = aligner.align(target_aa[target_order[position]], ref_aa[ref_order[position]])
        # Collect info needed - start and end positions of alignments ex. [0_old, 101]
        target_chain_info = [temp_align[0].path[0][0], temp_align[0].path[1][0]]
        ref_chain_info = [temp_align[0].path[0][1], temp_align[0].path[1][1]]
        # Collect XYZ cords. for target
        last_pos = 0
        count = -1
        for atom in target_atoms[target_order[position]]:
            if atom["atom_id"] == "CA" or not ca:  # Check if not carbon alpha
                if atom["atom_id"][0] != "H":  # Skip hydrogen atoms added while docking
                    if count == -1:
                        last_pos = atom["comp_num"]
                        count = 0
                    if atom["comp_num"] != last_pos:
                        count += 1
                    if count in range(target_chain_info[0], target_chain_info[1]):
                        rmsd_array["target"].append([atom["X"], atom["Y"], atom["Z"]])
        # Collect XYZ cords. for ref
        count = -1
        for atom in ref_atoms[ref_order[position]]:
            if atom["atom_id"] == "CA" or not ca:  # Check if not Carbon alpha
                if atom["atom_id"][0] != "H":  # Skip hydrogen atoms added while docking
                    if count == -1:
                        last_pos = atom["comp_num"]
                        count = 0
                    if atom["comp_num"] != last_pos:
                        count += 1
                    if count in range(ref_chain_info[0], ref_chain_info[1]):
                        rmsd_array["ref"].append([atom["X"], atom["Y"], atom["Z"]])
    # print(rmsd_array)
    # print(len(rmsd_array["target"]))
    # print(len(rmsd_array["ref"]))
    # Calculate RMSD
    sum_ = 0
    for position in range(len(rmsd_array["target"])):
        tar = rmsd_array["target"][position]
        ref = rmsd_array["ref"][position]
        sum_ += (pow(tar[0] - ref[0], 2) + pow(tar[1] - ref[1], 2) + pow(tar[2] - ref[2], 2))
    _rmsd = sqrt(sum_ / len(rmsd_array["target"]))
    logging.info(f"Calculated RMSD: {_rmsd}")
    return _rmsd


def get_center(pdb: str) -> list:
    """
    Calculates the mean of all atoms in the PDB.

    :param pdb: Path to the PDB to be read
    :return: A list consisting of three entries which correspond to the mean
             of all atom coordinates in the X, Y, Z dimensions, in that order
    """
    logging.info(f"Calculating center coordinates for '{pdb}'...")
    default_atoms = {}
    coords_dict = {'X': [], 'Y': [], 'Z': []}
    for chain in get_chains(pdb):
        atoms = get_atoms_on_chain(pdb, chain)
        if len(atoms) == 0:
            logging.error(f"Couldn't get atoms for chain {chain} in '{pdb}'. Proceeding anyway will likely cause the "
                          f"resultant PDB to be incomplete.")
            if not cm().get("permissive"):
                sys.exit(1)
        default_atoms[chain] = atoms
        for atom in default_atoms[chain]:
            for key in coords_dict.keys():
                coords_dict[key].append(atom[key])
    # Calculate center
    centroid = [round(statistics.mean(coords_dict['X']), 3), round(statistics.mean(coords_dict['Y']), 3),
                round(statistics.mean(coords_dict['Z']), 3)]
    logging.info(f"Centroid for '{pdb}' is {centroid}")
    return centroid


def center(pdb: str, rename: str = None) -> None:
    """
    Gather all atoms in the file and calculate the average XYZ coord to find the center of the structure. Write out
    centered PDB.

    :param pdb: Path to the PDB to be read
    :param rename : Optional name for new PDB created by method with centered structure
    """
    logging.info(f"Trying to find center of structure in '{pdb}'...")
    header = ""
    with open(pdb, "r+") as pdb_file:
        for line in pdb_file:
            if line[0:6] in ["HEADER", "REMARK", "SSBOND"]:
                header += line
            elif line[0:6] == "ATOM  ":  # We've moved past the header
                break
    atoms = []
    full_atom = []
    # Collect atom information from PDB
    for chain in get_chains(pdb):
        atoms = get_atoms_on_chain(pdb, chain)
        if len(atoms) == 0:
            logging.error(f"Couldn't get atoms for chain {chain} in '{pdb}'. Proceeding anyway will likely cause the "
                          f"resultant PDB to be incomplete or cause other unexpected behaviour.")
            if not cm().get("permissive"):
                sys.exit(1)
        chain_atoms = atoms
        for atom in chain_atoms:
            # Append atom_line with full atom information
            full_atom.append(atom)
            # Append atoms with just the coordinates of each atom
            atoms.append([atom['X'], atom['Y'], atom['Z']])
    # Convert to array
    atom_array = np.array(atoms)
    # Calculate avg for each axis
    current_avg = np.mean(atom_array, axis=0)
    # Center the coordinates
    new_array = []
    for set_pos in atom_array:
        x = set_pos - current_avg
        new_array.append(x)
    # Calculate PCA
    pca = PCA(n_components=3)
    pca.fit(new_array)
    transpose = np.transpose(pca.components_)
    # Determinate of components matrix - Corrects if determinate is negative
    determinate = np.linalg.det(transpose)
    if determinate < 0:
        for position in transpose:
            position[0] = position[0] * -1
    # Rotate x-axis
    # Always have N-terminus in positive coordinates (Fixes flips on x-axis)
    test_x = np.matmul(new_array[-1], transpose)
    if test_x[1] < 0:
        r = Rotation.from_euler('x', 180, degrees=True)
        transpose = r.apply(transpose)
    # Rotate y-axis
    # Always have first chain on left side | alpha on left side (Fixes flips on y-axis)
    # Determines by looking at last atom which should be on the Beta chain
    test_y = np.matmul(new_array[-1], transpose)
    if test_y[0] < 0:
        r = Rotation.from_euler('y', 180, degrees=True)
        transpose = r.apply(transpose)
    # Multiply by EV
    new_cords = np.matmul(new_array, transpose)
    # Replace XYZ coordinates
    axis = ['X', 'Y', 'Z']
    for num in range(0, len(full_atom)):
        for position in range(0, len(axis)):
            full_atom[num][axis[position]] = new_cords[num][position]
    # Write to new file
    new_name = pdb[:-4] + "_center.pdb"
    if rename:
        new_name = rename
    with open(new_name, "w") as f:
        f.write(header)
        f.write(rebuild_atom_line(full_atom))
        f.write("TER\nEND\n")
        logging.info(f"Wrote centered pdb to {new_name}")


# TODO: Have this keep HEADER, REMARK, and SSBOND records as well.
#  Should this be taken from pdb_1 or pdb_2? Configurable?
def join(pdb_1: str, pdb_2: str, rename: str = None) -> None:
    f"""
    Joins together two PDB files by appending first PDBs atoms to second PDBs atoms

    :param pdb_1 : Location and name of PDB 1
    :param pdb_2 : Location and name of PDB 2
    :param rename : Name of new file (default "{pdb_1}_{pdb_2}_concat.pdb"
    """
    logging.info(f"Trying to join '{pdb_1}' and '{pdb_2}'...")
    atoms_lines = []
    pdbs = [pdb_1, pdb_2]
    for pdb in pdbs:
        with open(pdb, "r") as i:
            for line in i:
                if line[0:6] == "ATOM  " or line[0:6] == "TER   ":
                    atoms_lines.append(line)
    new_name = pdb_1[:-4] + "_" + pdb_2.split(os.sep)[-1][:-4] + "_concat.pdb"
    if rename:
        new_name = rename
    with open(new_name, "w") as o:
        for line in atoms_lines:
            o.write(line)
        o.write("END\n")
        logging.info(f"Wrote joined PDB to '{new_name}'!")


# TODO: Have this keep HEADER, REMARK, and SSBOND records as well.
def reorder_chains(pdb: str, chain_order: str, rename: str = None) -> None:
    """
    Update the chain order. Must send in a list with identical number of chains

    :param pdb: The path to the PDB to be read
    :param chain_order: Order in which you want the chains to be in the file
    :param rename: Optional name of the newly written PDB. Default is "{pdb}_reordered.pdb"
    """
    logging.info(f"Reordering chains in '{pdb}' with order '{chain_order}'...")
    chain_info = {}
    new_order = []
    for chain in get_chains(pdb):
        atoms = get_atoms_on_chain(pdb, chain)
        if len(atoms) == 0:
            logging.error(f"Couldn't get atoms for chain {chain} in '{pdb}'. Proceeding anyway will likely cause the "
                          f"resultant PDB to be incomplete.")
            if not cm().get("permissive"):
                sys.exit(1)
        chain_info[chain] = atoms
    for chain in list(chain_order):
        for atom in chain_info[chain]:
            new_order.append(atom)
    new_name = pdb[:-4] + "_reordered.pdb"
    if rename:
        new_name = rename
    with open(new_name, "w") as f:
        f.write(rebuild_atom_line(new_order))
        f.write("TER\nEND\n")
        logging.info(f"Wrote reordered PDB to '{new_name}'!")


# TODO: Have this keep HEADER, REMARK, and SSBOND records as well.
def update_chain_id(pdb: str, id_mapping: dict, rename: str = None) -> None:
    """
    Update the labels based on submitted chain dictionary
    Ex of dictionary: {'A':'D', 'B':'E'}  A gets replaced with D and B gets replaced with E

    :param pdb: Path to the PDB to be read
    :param id_mapping: Dictionary of chain ids that need to be adjusted. Every chain id of _key_ gets replaced with _val_
    :param rename: Optional new name for the PDB that is written. Default is '{pdb}_upd_chains.pdb'
    """
    logging.info(f"Updating chain ids in '{pdb}' using mapping {id_mapping}...")
    chain_info = {}
    new_order = []
    for chain in get_chains(pdb):
        atoms = get_atoms_on_chain(pdb, chain)
        if len(atoms) == 0:
            logging.error(f"Couldn't get atoms for chain {chain} in '{pdb}'. Proceeding anyway will likely cause the "
                          f"resultant PDB to be incomplete.")
            if not cm().get("permissive"):
                sys.exit(1)
        chain_info[chain] = atoms
    for chain in chain_info:
        for atom in chain_info[chain]:
            atom['chain_id'] = id_mapping[chain]
            new_order.append(atom)
    new_name = pdb[:-4] + "_upd_chains.pdb"
    if rename:
        new_name = rename
    with open(new_name, "w") as f:
        f.write(rebuild_atom_line(new_order))
        f.write("TER\nEND\n")
        logging.info(f"Wrote PDB with updated chain ids to '{new_name}'")


def match_number(pdb: str, upd_chains: str, pdb_ref: str, custom: str = None, rename: str = None) -> Path:
    """
    Align chains and for amino acids that overlap, renumber them to the reference chains
    numbering. Writes an updated PDB to disk.

    :param pdb: The PDB whose chains shall be aligned
    :param upd_chains: Which chains in the base PDB shall be aligned, should be uppercase chain ids without separators
    :param pdb_ref: The reference PDB to align against
    :param custom: Optional comma separated str of numbering to apply to the shortest chain
    :param rename: Optional name for the updated PDB. Default is '{pdb}_orinum.pdb'
    :return: A Path object to the written out PDB
    """

    def find_gap(seq: str) -> int:
        """
        Counts the number of gap characters ('-') in seq from the start until the first non-gap character appears.

        :param seq: sequence string to be searched for non-gap character ("-")
        :return Index of first non-gap character
        """
        count = 0
        for letter in seq:
            if letter == "-":
                count += 1
            else:
                break
        return count

    ref_chains = get_chains(pdb_ref)
    ref_seqs = {}  # Dictionary of seq of ref_seqs chain id: seq
    target_chains = get_chains(pdb)
    target_seqs = {}  # Dictionary containing chain id: seq
    ref_aa = {}  # Dictionary containing aa numbering of ref.

    header = (
        f"EXPDTA     MODEL                  RENUMBERED\nREMARK    3                                 \nREMARK    3 "
        f"PROGRAM    : VIPER, PDBtool {TOOL_VER} \n")
    ssbond = ""
    ssbond_res_ids = {}  # Will become mapping for updated ids

    for chain in target_chains:
        ref_seqs[chain] = get_amino_acids_on_chain(pdb_ref, chain)
        target_seqs[chain] = get_amino_acids_on_chain(pdb, chain)
        if chain in ref_chains:
            ref_aa[chain] = collections.OrderedDict()  # Dictionary for num aa: aa of target.
    with open(pdb_ref, "r") as r:
        for line in r:
            if line[0:6] == 'ATOM  ':  # Only references atoms
                if line[section_chain_id] in ref_aa.keys():
                    ref_aa[line[section_chain_id]][
                        int(line[section_residue_number[0]:section_residue_number[1]])] = line[section_residue[0]:
                                                                                               section_residue[1]]
    aligns = {}  # alignments to determine start of chain number if they don't start at the same point
    aligner = Bio.Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -.5
    for chain in target_chains:
        aligns[chain] = aligner.align(target_seqs[chain], ref_seqs[chain])
        # Determine best alignment
        curr_best = None
        for a in aligns[chain]:
            if curr_best is None:
                curr_best = a
                continue
            if a.score > curr_best.score:
                curr_best = a
        # If there is a gap at the start of the target seq, removes those positions from ref_aa until it matches
        start_of_chain = find_gap(aligns[chain].sequences[0])
        while start_of_chain > 0:
            ref_aa[chain].popitem(last=False)
            start_of_chain -= 1
    shortest_num = 10000
    shortest_chain = ""
    # For custom numbering, determining shortest chain
    for chain in target_seqs:
        if len(target_seqs[chain]) < shortest_num:
            shortest_chain = chain
            shortest_num = len(target_seqs[chain])
    if custom is not None:
        custom_num = custom.split(",")
        ref_aa[shortest_chain] = {}
        for aa in target_seqs[shortest_chain]:
            ref_aa[shortest_chain][custom_num.pop(0)] = one_to_three(aa)
    body = ""
    with open(pdb, "r") as t:
        atom_count = 1  # No matter what, starts atom count at 1
        current_residue_count = -10000
        current_chain = ""
        for line in t:
            if line[0:6] == "HEADER" or line[0:6] == "REMARK":
                header += line
                continue
            if line[0:6] == "SSBOND":
                # Prepare SSBOND mapping, save current SSBOND values
                ssbond += line
                ssbond_res_ids[(int(line[section_ssbond_residue_id1[0]:section_ssbond_residue_id1[1]]),
                                line[section_ssbond_chain_id1])] = -1
                ssbond_res_ids[(int(line[section_ssbond_residue_id2[0]:section_ssbond_residue_id2[1]]),
                                line[section_ssbond_chain_id2])] = -1
                continue
            if line[0:6] == 'ATOM  ':
                pdb_chain_id = line[section_chain_id]
                # This is the case at the very beginning
                if current_residue_count < 0:
                    current_residue_count = int(line[section_residue_number[0]:section_residue_number[1]])
                    current_chain = pdb_chain_id  # Tells us the starting chain
                # Checks to see if we're on a new chain
                if current_chain != pdb_chain_id:
                    ref_aa.pop(current_chain)  # Removes previous chain from list
                    current_chain = pdb_chain_id
                    current_residue_count = int(line[section_residue_number[0]:section_residue_number[1]])
                # We're looking at a new residue and must therefore update our current_residue_count
                if current_residue_count != int(line[section_residue_number[0]:section_residue_number[1]]):
                    # Removes current (leftmost) reference residue
                    ref_aa[pdb_chain_id].pop([*ref_aa[pdb_chain_id].keys()][0])
                    current_residue_count = int(line[section_residue_number[0]:section_residue_number[1]])
                # Update SSBOND mapping to new num
                for k in ssbond_res_ids.keys():
                    if current_residue_count == k[0] and pdb_chain_id == k[1]:
                        # Save residue id of current (leftmost) reference residue as new id for SSBOND entry
                        ssbond_res_ids[k] = [*ref_aa[pdb_chain_id].keys()][0]
                # Save updated line up to the updated atom count
                temp_line = "" + line[0:6] + str(atom_count).rjust(5)
                if pdb_chain_id.upper() in upd_chains:  # Skips chains that we're not renumbering
                    temp_line += line[11:section_residue_number[0]]  # Replaces gap from num to residue num
                    # This is the reference residue id for the current residue in the PDB
                    next_num = [*ref_aa[pdb_chain_id].keys()][0]
                    temp_line += str(next_num).rjust(4) + line[section_residue_number[1]:]
                    body += temp_line + "\n"
                else:  # If we're skipping a chain we just write it but keep atom numbering
                    temp_line += line[11:]
                    body += temp_line + "\n"
                atom_count += 1
                continue
            else:
                body += line
                continue
    # Update SSBOND entries based on the discovered mapping
    for old_entry, new_entry in ssbond_res_ids.items():
        old = str(old_entry[1]) + " " + str(old_entry[0]).rjust(4)
        new = str(old_entry[1]) + " " + str(new_entry).rjust(4)
        ssbond = ssbond.replace(old, new, 1)
    # Write out updated PDB
    new_name = pdb[:-4] + "_orinum.pdb"
    if rename:
        new_name = rename
    with open(new_name, "w") as w:
        w.write(header)
        w.write(ssbond)
        w.write(body)
    return Path(new_name)


def get_centroid(pdb: Union[str, Path], residues: List[Union[REBprocessor.Node, int]], weighted: bool = False,
                 to_chain: str = False) -> Tuple[float, float, float]:
    """
    Returns the centroid of the list of residues as a tuple of X, Y, Z coordinates. Can be weighted according to the
    interaction energies. If weighting is desired, a list of REBprocessor.Node instances and to_chain must be supplied.
    Interaction strengths get normalized to [0.001, 1] and will be multiplied with the difference of the centroid of
    each residue and the global (unweighted) centroid. If the interaction is repulsive, the centroid will be nudged
    away from the residue centroid by this multiplied difference, if it is attractive, it will be nudged towards the
    residue centroid by this multiplied difference.

    :param pdb: The path to the PDB file to be read
    :param residues: A list of either residue ids or REBprocessor.Node instances
    :param weighted: Whether to calculate a weighted centroid
    :param to_chain: If calculating a weighted centroid, which interaction energies to which chain to consider
    :return:
    """
    if isinstance(residues, list) and len(residues) == 0:
        if cm().get("permissive"):
            logging.warning(f"Passed empty list to get_centroid()! Returning (0.0, 0.0, 0.0)...")
            return 0.0, 0.0, 0.0
        else:
            logging.error(f"Cannot pass empty list to get_centroid()! Aborting...")
            raise ValueError("Cannot pass empty list to get_centroid()!")
    if weighted and not to_chain:
        if cm().get("permissive"):
            logging.warning(f"Trying to calculate weighted centroid, but 'to_chain' has not been specified! "
                            f"Returning unweighted centroid instead...")
            weighted = False
        else:
            logging.error("Trying to calculate weighted centroid, but 'to_chain' has not been specified! Aborting...")
            raise ValueError("Trying to calculate weighted centroid, but 'to_chain' has not been specified! "
                             "Aborting...")
    if weighted and (isinstance(residues[0], REBprocessor.Node) and residues[0].chain == to_chain):
        logging.warning("Using strength to own chain, this may cause unexpected behavior!")
    posx = 0.0
    posy = 0.0
    posz = 0.0
    ids = []
    strengths = {}
    if isinstance(residues, list):
        for r in residues:
            if weighted and not isinstance(r, REBprocessor.Node):
                if cm().get("permissive"):
                    logging.warning(f"Entry in residue list ({r}) wasn't a REBprocessor.Node! Skipping this one...")
                else:
                    logging.error(f"Entry in residue list ({r}) wasn't a REBprocessor.Node! Aborting...")
                    raise ValueError("Not all entries in the residue list passed to get_centroid() were "
                                     "REBprocesser.Node instances! This is necessary for weighted calculation!")
            ids.append(int(r))
            if weighted:
                strengths[r.residue_id] = r.strength.get(to_chain, 0)
    else:
        ids.append(residues)
    if len(ids) == 0:
        if cm().get("permissive"):
            logging.warning(f"Couldn't determine any valid residue ids! Returning (0_old, 0_old, 0_old)...")
            return 0.0, 0.0, 0.0
        else:
            logging.error(f"Couldn't determine any valid residue ids! Aborting...")
            raise ValueError("Couldn't determine any valid residue ids!")
    residue_pos = {}
    with open(pdb, "r+") as pdb_in:
        for line in pdb_in:
            if line[0:6] == "ATOM  " and int(line[section_residue_number[0]:section_residue_number[1]]) in ids:
                if int(line[section_residue_number[0]:section_residue_number[1]]) not in residue_pos:
                    residue_pos[int(line[section_residue_number[0]:section_residue_number[1]])] = [0.0, 0.0, 0.0, 0]
                residue_pos[int(line[section_residue_number[0]:section_residue_number[1]])][3] += 1
                residue_pos[int(line[section_residue_number[0]:section_residue_number[1]])][0] += float(
                    line[section_x[0]:section_x[1]])
                residue_pos[int(line[section_residue_number[0]:section_residue_number[1]])][1] += float(
                    line[section_y[0]:section_y[1]])
                residue_pos[int(line[section_residue_number[0]:section_residue_number[1]])][2] += float(
                    line[section_z[0]:section_z[1]])
    if not weighted:
        for pos in residue_pos.values():
            posx += pos[0] / pos[3]
            posy += pos[1] / pos[3]
            posz += pos[2] / pos[3]
        return posx / len(residue_pos), posy / len(residue_pos), posz / len(residue_pos)
    else:
        # TODO: Just avg residue position * normed strength instead of nudging routine?
        #  i. e.: ( (0_old.1 * xyz_1 + 0_old.4 * xyz_2 + 0_old.2 * xyz_3 + ...) / num_res )
        # Normalize values
        min_strength_attractive = 10000.0
        min_strength_repulsive = 10000.0
        max_strength_attractive = 0.0
        max_strength_repulsive = 0.0
        for s in strengths.values():
            if s < 0:
                if -1 * s < min_strength_attractive:
                    min_strength_attractive = -1 * s
                if -1 * s > max_strength_attractive:
                    max_strength_attractive = -1 * s
            elif s > 0:
                if s < min_strength_repulsive:
                    min_strength_repulsive = s
                    # min_strength_repulsive = s / 3
                if s > max_strength_repulsive:
                    max_strength_repulsive = s
                    # max_strength_repulsive = s / 3

        def _norm(strength: float) -> Tuple[float, int]:
            if strength < 0:
                return max(0.001, (-1 * s - min_strength_attractive) /
                           (max_strength_attractive - min_strength_attractive)), 1
            elif strength > 0:
                return max(0.001, (s - min_strength_repulsive) / (max_strength_repulsive - min_strength_repulsive)), -1
            else:
                return 0.0001, 1  # epsilon

        for rid, s in strengths.items():
            strengths[rid] = _norm(s)

        centroid_x = 0.0
        centroid_y = 0.0
        centroid_z = 0.0

        for pos in residue_pos.values():
            centroid_x += pos[0] / pos[3]
            centroid_y += pos[1] / pos[3]
            centroid_z += pos[2] / pos[3]

        centroid_x /= len(residue_pos)
        centroid_y /= len(residue_pos)
        centroid_z /= len(residue_pos)

        x_weighted = centroid_x
        y_weighted = centroid_y
        z_weighted = centroid_z

        for rid, s in strengths.items():
            # Add weighted difference of res_pos <> centroid to start centroid pos
            if s[1] == 1:
                x_weighted += (residue_pos[rid][0] / residue_pos[rid][3] - centroid_x) * s[0]
                y_weighted += (residue_pos[rid][1] / residue_pos[rid][3] - centroid_y) * s[0]
                z_weighted += (residue_pos[rid][2] / residue_pos[rid][3] - centroid_z) * s[0]
            elif s[1] == -1:  # in case of repulsion, push weighted center into opposite direction
                x_weighted -= (residue_pos[rid][0] / residue_pos[rid][3] - centroid_x) * s[0]
                y_weighted -= (residue_pos[rid][1] / residue_pos[rid][3] - centroid_y) * s[0]
                z_weighted -= (residue_pos[rid][2] / residue_pos[rid][3] - centroid_z) * s[0]

        return x_weighted, y_weighted, z_weighted


def get_dist_centroid(pdb: str, first: List[Union[int, REBprocessor.Node]], second: List[Union[int, REBprocessor.Node]],
                      resolution_total: bool = True) -> Union[float, Tuple[float, int, int]]:
    """
    Returns the distance between the centroid of the first and second list of residues. If resolution_total is set to
    False, it will instead calculate the distance between the closest two residues of the two lists and return a tuple
    of (distance, resid_from_first, resid_from_second).

    :param pdb: The path to the PDB to be read
    :param first: List of residue ids
    :param second: List of residue ids
    :param resolution_total: Whether to consider the aggregate centroid of all residues in one list,
            or consider the centroid of each residue itself
    :return: The distance between the centroid of the first and second list of residues
    """
    if len(first) == 0 or len(second) == 0:
        if cm().get("permissive"):
            return 0.0
        else:
            raise ValueError("Can't pass empty list to avg distance function.")
    first_ids = []
    second_ids = []
    for entry in first:
        first_ids.append(int(entry))
    for entry in second:
        second_ids.append(int(entry))
    x_avg1, y_avg1, z_avg1, x_avg2, y_avg2, z_avg2 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    num1, num2 = 0, 0
    last_seen = -1
    centroid_dict1 = {}
    centroid_dict2 = {}
    with open(pdb, "r") as read:
        for line in read:
            if line[0:6] != "ATOM  ":
                continue
            else:
                residue_number = int(line[section_residue_number[0]:section_residue_number[1]])
                if residue_number != last_seen:
                    last_seen = residue_number
                if residue_number in first_ids:
                    if not resolution_total:
                        if c := centroid_dict1.get(str(residue_number), None):
                            c[0] += float(line[section_x[0]:section_x[1]])
                            c[1] += float(line[section_y[0]:section_y[1]])
                            c[2] += float(line[section_z[0]:section_z[1]])
                            c[3] += 1
                        else:
                            centroid_dict1[str(residue_number)] = [float(line[section_x[0]:section_x[1]]),
                                                                   float(line[section_y[0]:section_y[1]]),
                                                                   float(line[section_z[0]:section_z[1]]),
                                                                   1]
                    else:
                        num1 += 1
                        x_avg1 += float(line[section_x[0]:section_x[1]])
                        y_avg1 += float(line[section_y[0]:section_y[1]])
                        z_avg1 += float(line[section_z[0]:section_z[1]])
                elif residue_number in second_ids:
                    if not resolution_total:
                        if c := centroid_dict2.get(str(residue_number), None):
                            c[0] += float(line[section_x[0]:section_x[1]])
                            c[1] += float(line[section_y[0]:section_y[1]])
                            c[2] += float(line[section_z[0]:section_z[1]])
                            c[3] += 1
                        else:
                            centroid_dict2[str(residue_number)] = [float(line[section_x[0]:section_x[1]]),
                                                                   float(line[section_y[0]:section_y[1]]),
                                                                   float(line[section_z[0]:section_z[1]]),
                                                                   1]
                    else:
                        num2 += 1
                        x_avg2 += float(line[section_x[0]:section_x[1]])
                        y_avg2 += float(line[section_y[0]:section_y[1]])
                        z_avg2 += float(line[section_z[0]:section_z[1]])
    if resolution_total:
        return sqrt((x_avg1 / num1 - x_avg2 / num2) ** 2 +
                    (y_avg1 / num1 - y_avg2 / num2) ** 2 +
                    (z_avg1 / num1 - z_avg2 / num2) ** 2)
    else:
        # TODO: Is this faster with a k-d tree?
        dists = []
        for res1, centroid1 in centroid_dict1.items():
            for res2, centroid2 in centroid_dict2.items():
                dists.append([sqrt((centroid1[0] / centroid1[3] - centroid2[0] / centroid2[3]) ** 2 +
                                   (centroid1[1] / centroid1[3] - centroid2[1] / centroid2[3]) ** 2 +
                                   (centroid1[2] / centroid1[3] - centroid2[2] / centroid2[3]) ** 2),
                              res1, res2])
        closest_centroids = min(dists, key=lambda d: d[0])
        return closest_centroids[0], int(closest_centroids[1]), int(closest_centroids[2])


_tree_cache = {}


def get_dist_closest_atom(pdb: str, first: List[Union[int, REBprocessor.Node]],
                          second: List[Union[int, REBprocessor.Node]]) -> Tuple[float, int, int]:
    """
    Returns the distance between the closest two atoms of the closest two residues in first and second, as well as which
    residues are involved from each list.

    :param pdb: The path to the PDB to be read
    :param first: A list of residue ids
    :param second: A list of residue ids
    :return: (closest distance, residue 1, residue 2)
    """
    if len(first) == 0 or len(second) == 0:
        if cm().get("permissive"):
            return 0.0, 0, 0
        else:
            raise ValueError("Can't pass empty list to avg distance function.")
    first_ids = []
    second_ids = []
    for entry in first:
        first_ids.append(int(entry))
    for entry in second:
        second_ids.append(int(entry))

    use_tree = _tree_cache.get(frozenset(first_ids), None)
    if len(_tree_cache) > 20:  # Free up memory
        _tree_cache.clear()

    points_first = []
    points_second = []
    with open(pdb, "r") as read:
        for line in read:
            if line[0:6] != "ATOM  ":
                continue
            else:
                residue_number = int(line[section_residue_number[0]:section_residue_number[1]])
                if residue_number in first_ids:
                    points_first.append([float(line[section_x[0]:section_x[1]]),
                                         float(line[section_y[0]:section_y[1]]),
                                         float(line[section_z[0]:section_z[1]]),
                                         residue_number])
                elif residue_number in second_ids:
                    points_second.append([float(line[section_x[0]:section_x[1]]),
                                          float(line[section_y[0]:section_y[1]]),
                                          float(line[section_z[0]:section_z[1]]),
                                          residue_number])
    if use_tree is None:
        use_tree = KDTree([p[0:-1] for p in points_first])
        _tree_cache[frozenset(first_ids)] = use_tree

    distances, neighbor_indices = use_tree.query([p[0:-1] for p in points_second], k=1)
    closest = distances.min()
    distances, neighbor_indices = distances.tolist(), neighbor_indices.tolist()

    idx_second = distances.index(closest)
    idx_first = neighbor_indices[idx_second][0]

    return closest, points_first[idx_first][-1], points_second[idx_second][-1]

