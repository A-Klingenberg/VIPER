# This code has been adapted from previous work by Austin Seamann, Dario Ghersi, and Ryan Ehrlich.
# Please refer to https://github.com/Aseamann/ACEdecoy

import logging
import os
import sys
import statistics
from math import sqrt
from typing import Optional, List

import Bio
import numpy as np
from Bio import Align, pairwise2
from scipy.spatial.transform import Rotation
from sklearn.decomposition import PCA

from VIPER import configmanager as cm

TOOL_VER = "2.0"

# Magic numbers to extract data from PDB format.
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


def get_chains(pdb: str) -> Optional[list]:
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
        if cm.get("permissive"):
            return None
        else:
            sys.exit(1)
    return chains


def get_amino_acids_on_chain(pdb: str, chain: str) -> str:
    """
    Returns the string of amino acids in a specific chain as a string in single letter notation.

    :param pdb: Path to PDB file to read
    :param chain: Which chain to read
    :return: String of amino acids in single letter formatting, or None, if chain can't be found in PDB and VIPER is
        running in permissive mode
    """
    logging.debug(f"Trying to get amino acid sequence of chain '{chain}' in '{pdb}'...")
    output = ''
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
                            output += three_to_one(line[section_residue[0]:section_residue[1]])
                            count += 1
                    elif count < int(line[section_residue_number[0]:section_residue_number[1]]):
                        count = int(line[section_residue_number[0]:section_residue_number[1]])
    logging.debug(f"Got following sequence: '{output}'")
    if len(output) == 0:
        logging.error(f"Couldn't find any amino acid in chain '{chain}' in '{pdb}'!")
        if not cm.get("permissive"):
            sys.exit(1)
    return output


def three_to_one(three: str) -> str:
    """
    Converts three letter amino acid abbreviation to single letter abbreviation.

    :param three: Three letter amino acid abbreviation
    :return: Single letter amino acid abbreviation, or ' ' if input abbreviation is not recognized and VIPER is running
        in permissive mode.
    """
    translate = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'ASX': 'B', 'CYS': 'C', 'GLU': 'E',
        'GLN': 'Q', 'GLX': 'Z', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y',
        'VAL': 'V'
    }
    if three.upper() in translate:
        return translate[three.upper()]
    else:
        logging.log(logging.WARN, f"Tried to convert '{three}' to single letter amino acid abbreviation but failed!")
        return " "


def one_to_three(one: str) -> str:
    """
    Converts three letter amino acid to single letter abbreviation.

    :param one: One letter amino acid abbreviation
    :return: Three letter amino acid abbreviation. If the one letter amino acid abbreviation is not recognized, an VIPER
        is run in permissive mode, return '   ', otherwise stop program.
    """
    translate = {
        'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'B': 'ASX', 'C': 'CYS', 'E': 'GLU',
        'Q': 'GLN', 'Z': 'GLX', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
        'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR',
        'V': 'VAL'
    }
    if one.upper() in translate:
        return translate[one.upper()]
    else:
        logging.log(logging.WARN, f"Tried to convert '{one}' to three letter amino acid abbreviation but failed!")
        return "   "


def first_atom_on_chain(pdb: str, chain: str) -> Optional[dict]:
    """
    Returns an atom dictionary for the first atom of a chain.

    :param pdb: Path to PDB to be read
    :param chain: Which chain to read
    :return: A dict of properties  of the first atom in the specified chain, or None, if no first atom could be found and
        VIPER is running in permissive mode
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
        if cm.get("permissive"):
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
        logging.error(f"Couldn't locate atom number {atom_num} from '{pdb}'!")
        if cm.get("permissive"):
            return None
        else:
            sys.exit(1)


def get_atoms_on_chain(pdb: str, chain: str) -> list:
    """
    Collect the atoms from an inputted chain. Provides all values that are in the specified PDB file.

    :param pdb: Path to PDB to be read
    :param chain: Which chain to use
    :return: List of atoms based on the chain id that was submitted. This list may be empty if no atoms could be found
        and VIPER is running permissive mode
    """
    logging.debug(f"Trying to get atoms on chain '{chain}' in '{pdb}'...")
    atoms = []
    with open(pdb, 'r') as file:
        for line in file:
            if line[0:6] == 'ATOM  ':
                if line[section_chain_id] == chain and len(line) >= 76:
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
                    atoms.append(atom)
    if len(atoms) == 0:
        logging.error(f"Tried to get atoms on chain '{chain}' from '{pdb}' but failed!")
        if not cm.get("permissive"):
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
    :return: Distance between two atoms in angstroms, or None, if the distance could not be computed and VIPER is running
        in permissive mode
    """
    logging.debug(f"Trying to find euclidean distance between atom {atom_num_1} and {atom_num_2} in '{pdb}'...")
    atom_1 = get_atom(pdb, atom_num_1)
    atom_2 = get_atom(pdb, atom_num_2)
    if atom_1 is None or atom_2 is None:
        logging.error(
            f"Tried to compute euclidean distance between atom {atom_num_1} and {atom_num_2} in '{pdb}' but failed!")
        if cm.get("permissive"):
            return None
        else:
            sys.exit(1)
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
    logging.debug(f"Rebuilt atom(s): {output}")
    return output


def renumber_ascending(pdb: str, rename: str = None) -> None:
    """
    Write out a PDB file in the same directory as the input PDB with atom number and sequence number in ascending order.
    There is a limit on how large atom (serial) ids and residue ids may be inherent in the PDB format (99,999 and 9,999).
    If the proteins are too large, the count may exceed this number. Should VIPER be running in permissive mode, the ids
    will be taken modulo 100,000 and 10,000, respectively.
    This may cause errors later on, since ids will be duplicated!

    :param pdb: Path to PDB to be read
    :param rename: Optional name for PDB file being created, otherwise use "_renumbered" appended to original name
    """
    logging.info(f"Renumbering {pdb} by counting up...")

    out_name = pdb[:-4] + "_renumbered.pdb"
    if rename:
        out_name = rename

    running_atom_count = 1  # Keeps track of atom number
    previous_residue_number = -10000
    running_residue_count = 0  # Keeps track of residue number
    chains = []  # Chains already seen
    expdta_written = False
    with open(pdb, "r") as i, open(out_name, "w+") as o:
        for line in i:
            if line[0:6] == 'HEADER':
                o.write(line)
            if line[0:6] == 'ATOM  ':  # Only copy over atoms
                if not expdta_written:
                    o.write('EXPDTA    DOCKING MODEL           RENUMBERED\n')
                    o.write('REMARK    3                                 \n')
                    o.write(f'REMARK    3 PROGRAM    : VIPER, PDBtool {TOOL_VER} \n')
                    expdta_written = True
                if line[16] != 'B' and line[26] == ' ':  # Don't allow secondary atoms
                    current_atom_number = line[section_atom_number[0]:section_atom_number[1]]
                    current_residue_number = int(line[section_residue_number[0]:section_residue_number[1]])
                    if len(chains) == 0:
                        chains.append(line[section_chain_id])
                    elif line[section_chain_id] != chains[-1]:  # Encountered new chain
                        chains.append(line[section_chain_id])
                        o.write('TER\n')
                    if current_residue_number != previous_residue_number:  # Encountered next residue
                        previous_residue_number = int(line[section_residue_number[0]:section_residue_number[1]])
                        running_residue_count += 1
                    # Update atom count and residue count
                    # Because there are only 5 characters for the atom serial number and 4 for the residue, we need
                    # to restrict the number space to < 100,000 and 10,000, respectively
                    if running_atom_count > 99999 or running_residue_count > 9999:
                        if cm.get("permissive"):
                            logging.log(logging.WARN, f"While renumbering {pdb} atom or residue count exceeded legal "
                                                      f"limits. Counts will be taken modulo 100,000 and 10,000, "
                                                      f"respectively. Proceeding may cause unexpected program "
                                                      f"behaviour.")
                            line = line[:6] + str(running_atom_count % 100000).rjust(5) + line[11:22] + str(
                                running_residue_count).rjust(4) + line[26:]
                        else:
                            logging.log(logging.ERROR, f"While renumbering {pdb} atom or residue count exceeded legal "
                                                       f"limits (100000 or 10000, respectively. Aborting! You might "
                                                       f"want to try to reduce the proteins to a cut out version "
                                                       f"focusing on the binding site.")
                            sys.exit(1)
                    running_atom_count += 1
                    o.write(line)
        o.write('TER\nEND\n')


def remove_chain(pdb: str, chain_id: List[str], rename: str = None) -> None:
    """
    Removes chains with a certain ID and save the results in a PDB.

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
            if line[0:6] == 'ATOM  ' or line[0:6] == 'TER   ':
                if line[section_chain_id].upper() not in ids:
                    o.write(line)
        logging.info(f"PDB with chains {ids} removed saved to {rename}!")


# FIXME: Produces faulty PDB? When written out, there is data following TER in the same line
def superimpose(pdb: str, ref_pdb: str, target_order: str, ref_order: str, rename: bool = None) -> float:
    """
    Superimpose two PDBs.
    Read in target and reference structure, superimpose target to reference, save superimposed target structure.

    :param pdb: Location of target PDB
    :param ref_pdb: Location of reference PDB
    :param target_order: Order of chains to compare in superimposing structures for target
    :param ref_order: Order of chains to compare in superimposing structures for reference
    :param rename: Optional name in for resulting superimposed positioning of target structure
    :return RMSD value for all-atom RMSD
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.gap_score = -100.00
    aligner.match_score = 2.0
    aligner.mismatch_score = 0.0

    logging.info(f"Trying to superimpose {pdb} onto {ref_pdb} with chain order {target_order} (target) and "
                 f"{ref_order} (reference).")

    # Collect target pdb seq
    target_seq = {}  # Dictionary chain_id:seq
    for chain in get_chains(pdb):
        residues = get_amino_acids_on_chain(pdb, chain)
        if residues is None:
            logging.error(f"Couldn't find any residues for chain {chain} in '{pdb}'!")
            if not cm.get("permissive"):
                logging.error(f"Aborting now...")
                sys.exit(1)
        target_seq[chain] = residues

    # Collect reference pdb seq
    ref_seq = {}  # Dictionary chain_id:seq
    for chain in get_chains(ref_pdb):
        residues = get_amino_acids_on_chain(ref_pdb, chain)
        if residues is None:
            logging.error(f"Couldn't find any residues for chain {chain} in '{pdb}'!")
            if not cm.get("permissive"):
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
        logging.log(logging.INFO, f"Target Chain: {target_pos[chain_pos]}")
        logging.log(logging.INFO, f"Reference Chain: {ref_pos[chain_pos]}")
        logging.log(logging.INFO, f"{temp_align[0]}")

        start_pos['target'][target_pos[chain_pos]] = [temp_align[0].path[0][0], temp_align[0].path[1][0]]
        start_pos['reference'][ref_pos[chain_pos]] = [temp_align[0].path[0][1], temp_align[0].path[1][1]]

    # Initialize parser
    if cm.get("verbose"):
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
                    if first_in_chain:  # If chain doesn't start count at position 0
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
    new_name = pdb[:-4] + "_aligned.pdb"
    if rename:
        new_name = rename
    io.save(new_name)
    return super_imposer.rms


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
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.gap_score = -100.00
    aligner.match_score = 2.0
    aligner.mismatch_score = 0.0
    # Loop through each paired chain
    for position in range(len(target_order)):
        # run alignment
        temp_align = aligner.align(target_aa[target_order[position]], ref_aa[ref_order[position]])
        # Collect info needed - start and end positions of alignments ex. [0, 101]
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
            if not cm.get("permissive"):
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
    atoms = []
    full_atom = []
    # Collect atom information from PDB
    for chain in get_chains(pdb):
        atoms = get_atoms_on_chain(pdb, chain)
        if len(atoms) == 0:
            logging.error(f"Couldn't get atoms for chain {chain} in '{pdb}'. Proceeding anyway will likely cause the "
                          f"resultant PDB to be incomplete or other unexpected behaviour.")
            if not cm.get("permissive"):
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
        f.write(rebuild_atom_line(full_atom))
        logging.info(f"Wrote centered pdb to {new_name}")


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
        logging.info(f"Wrote joined PDB to '{new_name}'!")


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
            if not cm.get("permissive"):
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
        logging.info(f"Wrote reordered PDB to '{new_name}'!")


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
            if not cm.get("permissive"):
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
        logging.info(f"Wrote PDB with updated chain ids to '{new_name}'")


# FIXME: Also renumber SSBOND header entries
def match_number(pdb: str, upd_chains: str, pdb_ref: str, custom: str, rename: str = None) -> None:
    """
    Align chains and for amino acids that overlap, renumber them to the reference chains
    numbering. Update them to PDB.

    Parameters
    __________
    custom : str
        comma sep. str of numbering for shorest chain if provided by user
    """

    def find_gap(seq: str) -> int:
        """
        Counts how may gaps until start of seq and reports position

        Parameters
        __________
        seq : str
            sequence string to be search for gap character ("-")

        Returns
        _______
        count : int
            position in which first gap character appears
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

    for chain in target_chains:
        ref_seqs[chain] = get_amino_acids_on_chain(pdb_ref, chain)
        target_seqs[chain] = get_amino_acids_on_chain(pdb, chain)
        if chain in ref_chains:
            ref_aa[chain] = {}  # Dictionary for num aa: aa of target.
    with open(pdb_ref, "r") as r:
        for line in r:
            if line[0:6] == 'ATOM  ':  # Only references atoms
                if line[section_chain_id] in ref_aa.keys():
                    ref_aa[line[section_chain_id]][
                        int(line[section_residue_number[0]:section_residue_number[1]])] = line[section_residue[0]:
                                                                                               section_residue[1]]
    aligns = {}  # alignments to determine start of chain number if they don't start at the same point
    for chain in target_chains:
        aligns[chain] = pairwise2.align.globalms(target_seqs[chain], ref_seqs[chain], 2, -1, -2, -.5,
                                                 penalize_end_gaps=(False, False), one_alignment_only=True)
        # If there is a gap at the start of the target seq, removes those positions from ref_aa until it matches
        start_of_chain = find_gap(aligns[chain][0][0])
        while start_of_chain != 0:
            ref_aa[chain].pop(ref_aa[chain].key()[0])
            start_of_chain -= 1
    shortest_num = 10000
    shortest_chain = ""
    # For custom numbering, determining shortest chain
    for chain in target_seqs:
        if len(target_seqs[chain]) < shortest_num:
            shortest_chain = chain
            shortest_num = len(target_seqs[chain])
    if custom != "...":
        custom_num = custom.split(",")
        ref_aa[shortest_chain] = {}
        for aa in target_seqs[shortest_chain]:
            ref_aa[shortest_chain][custom_num.pop(0)] = one_to_three(aa)
    with open(pdb, "r") as t:
        new_name = pdb[:-4] + "_orinum.pdb"
        if rename:
            new_name = rename
        with open(new_name, "w") as w:
            atom_count = 1  # No matter what, starts atom count at 1
            aa_num = -10000
            current_chain = ""
            for line in t:
                if line[0:6] == 'ATOM  ':
                    if aa_num < 0:  # Checks to see if we're just starting
                        aa_num = int(line[section_residue_number[0]: section_residue_number[1]])  # Tell me current amino acid number
                        current_chain = line[section_chain_id]  # Tells us the starting chain
                    if current_chain != line[section_chain_id]:  # Checks to see if we're on a new chain
                        ref_aa.pop(current_chain)  # Removes previous chain from list
                        current_chain = line[section_chain_id]  # Updates current chain
                        aa_num = int(
                            line[section_residue_number[0]: section_residue_number[1]])  # Updates current aa count
                    if aa_num != int(line[section_residue_number[0]: section_residue_number[1]]):  # Tells us if we're on a new aa
                        ref_aa[line[section_chain_id]].pop(
                            list(ref_aa[line[section_chain_id]])[0])  # Removes last counted AA
                        aa_num = int(line[section_residue_number[0]: section_residue_number[1]])
                    temp_line = ""
                    temp_line += line[0:6]  # Adds header
                    temp_line += str(atom_count).rjust(5)  # Adds atom count
                    atom_count += 1
                    if line[section_chain_id].upper() in upd_chains:  # Skips chains that we're not renumbering
                        temp_line += line[11:section_residue_number[0]]  # Replaces gap from num to AA num
                        next_num = list(ref_aa[line[section_chain_id]])[0]
                        temp_line += str(next_num).rjust(4)
                        temp_line += line[section_residue_number[1]:]
                        w.write(temp_line)
                    else:  # If we're skipping a chain we just write it but we keep atom numbering
                        temp_line += line[11:]
                        w.write(temp_line)
                else:
                    w.write(line)
