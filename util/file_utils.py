import os
import zipfile
from pathlib import Path
from typing import Union, List


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


def gather_files(directory: Union[str, Path], type: str = "pdb", recursive: bool = False) -> List:
    """
    Returns a list of paths to all PDB files in a given directory (and its subdirectories if recursive is true).

    :param directory: Which directory to search for PDB files
    :param type: Which file type to gather. Default is 'pdb'
    :param recursive: Whether to also search subdirectories
    :return: A list of paths to PDB files
    """
    directory = Path(directory)
    if recursive:
        paths = []
        for file in directory.glob(f"**/*.{type}"):
            paths.append(file)
    else:
        paths = []
        for file in directory.glob(f"*.{type}"):
            paths.append(file)
    return paths
