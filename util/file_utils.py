import copy
import logging
import os
import uuid
import zipfile
from pathlib import Path
from typing import Union, List

from ConfigManager import ConfigManager

cm = ConfigManager.get_instance


def compress_directory(directory: Union[str, Path], archive_name: Union[str, Path],
                       compression_type=zipfile.ZIP_DEFLATED, exclude: Union[str, List] = "",
                       compresslevel: int = 9, delete_uncompressed: bool = True) -> None:
    """
    This method compresses all files in a directory and its subdirectories into a compressed archive.
    Can optionally exclude a list of files, and delete each original file that was compressed.
    Important note: Files in the directory with the same name as archive_name will automatically be excluded!

    :param directory: Which directory to search for files
    :param archive_name: Which name to give the compressed archive (full path to file included)
    :param compression_type: Which compression type in the zipfile module to use. Standard: DEFLATED
    :param exclude: A path or list of paths to files which shall not be included
    :param compresslevel: Which compression level to use. Standard: 9, best compression, worst speed
    :param delete_uncompressed: Whether to delete uncompressed files after they've been added to the archive
    :return: None
    """
    use_path = directory
    if not isinstance(directory, Path):
        use_path = Path(directory)
    if not isinstance(exclude, list):
        exclude = [exclude]
    # Make sure to not include the archive within itself, leading to infinite recursion.
    exclude.append(Path(archive_name).name)

    with zipfile.ZipFile(archive_name, "w", compression=compression_type, allowZip64=True,
                         compresslevel=compresslevel) as archive_out:
        for file in use_path.glob("**/*"):
            if file.name in exclude:
                continue
            archive_out.write(file, arcname=file.relative_to(use_path))
            if delete_uncompressed and os.path.isfile(file):
                os.remove(file)


def make_pdb_ensemble_list(directory: Union[str, Path], out_path: Union[str, Path]) -> str:
    """
    Writes a list of paths to all PDB files in a directory and its subdirectories into a file, one path per line.

    :param directory: Which directory to search for PDB files
    :param out_path: Into which file to write the paths to the PDB files
    :return: The path to the ensemble file (normalized 'out_path')
    """
    if isinstance(directory, str):
        directory = Path(directory)
    out = os.path.normpath(out_path)
    with open(out, "w") as ensemble_out:
        for file in directory.glob("**/*.pdb"):
            ensemble_out.write(os.path.normpath(file) + "\n")
    return out


def gather_files(directory: Union[str, Path], filetype: str = "pdb", recursive: bool = False) -> List:
    """
    Returns a list of paths to all PDB files in a given directory (and its subdirectories if recursive is true).

    :param directory: Which directory to search for PDB files
    :param filetype: Which file type to gather. Default is 'pdb'
    :param recursive: Whether to also search subdirectories
    :return: A list of paths to PDB files
    """
    directory = Path(directory)
    if recursive:
        paths = []
        for file in directory.glob(f"**/*.{filetype}"):
            paths.append(file)
    else:
        paths = []
        for file in directory.glob(f"*.{filetype}"):
            paths.append(file)
    return paths


def make_file(path: Union[str, Path, List[Union[str, Path]]], content: str, dry_run: bool = False) -> Path:
    """
    Writes a file with content 'content' to the results path, creating any subdirectories as needed.
    The path will be based off the results path, so do not use absolute paths!
    Instead, use a list with directory names and end it with the filename.

    :param path: Where to write the file to
    :param content: What to write to the file
    :param dry_run: Whether or not to actually write to the file - True means that the path will be validated and any
        necessary subdirectories created, but the file _won't_ actually be written out!
    :return: A Path object to the written out file
    """
    p = copy.deepcopy(path)
    if isinstance(path, str) or isinstance(path, Path):
        p = [copy.deepcopy(Path(path))]
    for index, elem in enumerate(p):
        if isinstance(elem, Path) and elem.is_absolute():
            logging.warning(f"You can't pass an absolute path ({elem}) to make_file()!")
            raise ValueError(f"You can't pass an absolute path ({elem}) to make_file()!")
        p[index] = Path(elem)
    fname = ""
    if p[-1].is_dir():
        if cm().get("permissive"):
            fname = str(uuid.uuid4())
            logging.warning(f"The passed path '{p}' is a directory! Using '{fname}' as filename.")
        else:
            logging.error(f"The path '{p}' leads to a directory! Can only write to a file!")
            raise ValueError(f"The path '{p}' leads to a directory! Can only write to a file!")
    use_path = os.path.normpath(cm().get("results_path"))
    for elem in p:
        use_path = os.path.join(use_path, elem)
    use_path = Path(os.path.normpath(use_path))
    os.makedirs(use_path.parents[0], exist_ok=True)
    if not dry_run:
        with open(use_path, "w") as f:
            f.write(content)
            logging.info(f"Wrote out file '{use_path}'")
    else:
        logging.info(f"Did dry run of make_file() for '{use_path}'")
    return use_path
