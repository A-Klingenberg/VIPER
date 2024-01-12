import logging
import re
import time

import requests
from requests_toolbelt.multipart.encoder import MultipartEncoder


from ConfigManager import ConfigManager

cm = ConfigManager.get_instance


class PEPstrMODWrapper:
    ENDPOINT = "https://webs.iiitd.edu.in/raghava/pepstrmod/nat_ss_upload.php"
    ENDPOINT_REDUNDANCY = "http://osddlinux.osdd.net/raghava/pepstrmod/nat_ss_upload.php"
    _use_redundancy = False

    @staticmethod
    def submit_peptide(sequence: str, options: dict = None) -> str:
        """
        Submits a peptide sequence to the PEPstrMOD service, optionally overriding the settings in the config.json,
        and returns the PDB as a str. This might take a few minutes.

        :param sequence: The peptide sequence in single letter format
        :param options: Any optional override settings
        :return: The PDB as a str
        """
        # TODO: Find out hardcoded values for all environment options
        form_data = {
            "MAX_FILE_SIZE": str(1000000),
            "seq": sequence,
            "email": cm().get("pepstrmod_config.email"),
            "time": cm().get("pepstrmod_config.simulation_time"),
            "env": cm().get("pepstrmod_config.environment"),
            "topol": "yes" if cm().get("pepstrmod_config.download_topology") else "no",
            "cluster": "yes" if cm().get("pepstrmod_config.cluster_analysis") else "no",
            "traj": "yes" if cm().get("pepstrmod_config.download_trajectory") else "no",
            "graph": "yes" if cm().get("pepstrmod_config.do_energy_rms_graph") else "no",
            "send": "Submit sequence for prediction",
        }

        if options:
            for k, v in options.items():
                form_data[k] = v

        # TODO: Do this by hand, so that the dependency can be dropped
        multipart_data = MultipartEncoder(
            fields={
                "MAX_FILE_SIZE": form_data["MAX_FILE_SIZE"],
                "seq": form_data["seq"],
                "email": form_data["email"],
                "time": form_data["time"],
                "env": form_data["env"],
                "topol": form_data["topol"],
                "cluster": form_data["cluster"],
                "traj": form_data["traj"],
                "graph": form_data["graph"],
                "send": form_data["send"],
            }
        )

        response = requests.post(url=PEPstrMODWrapper.ENDPOINT, data=multipart_data,
                                 headers={'Content-Type': multipart_data.content_type}, timeout=60)

        run_id = None
        # TODO: Think of better regex?
        pattern = re.compile(
            r'(https?://(www\.)?webs\.iiitd\.edu\.in/raghava/pepstrmod/results/index\.php\?ran=(?P<runid>\d+))')
        if matched_url := re.search(pattern, response.text):
            run_id = matched_url.group("runid")
            logging.debug(f"Identified PEPstrMOD run id {run_id}")
        else:
            raise ValueError("Could not identify where the PEPstrMOD results (run id) will be posted. The results "
                             "endpoint has likely changed, please create an issue for this on the GitHub repository.")

        pdb_url = f"https://webs.iiitd.edu.in/dis.php?file=/tmp/pepstrmod/pep_{run_id}/pepstr{run_id}.pdb"
        logging.debug(f"Using PEPstrMOD endpoint {pdb_url}")

        # TODO: Find better marker for the results being ready?
        results_page = requests.get(pdb_url, timeout=300)
        interval = 60
        wait_limit = 30
        wait_count = 0
        while "HEADER PEPstrMOD" not in results_page.text:
            if wait_count > wait_limit:
                logging.debug(f"Exceeded wait limit of {wait_limit*interval} seconds for official endpoint.")
                PEPstrMODWrapper._use_redundancy = True
                break
            # Wait for results to be ready
            wait_count += 1
            time.sleep(interval)
            results_page = requests.get(pdb_url,
                                        headers={"Cache-Control": "no-cache", "Pragma": "no-cache",
                                                 "Host": "webs.iiitd.edu.in", "Referer": results_page.url},
                                        allow_redirects=True, stream=False, timeout=300)

        # official timed out, use redundancy
        if PEPstrMODWrapper._use_redundancy:
            logging.warning("Official PEPstrMOD service timed out, trying to use backup!")
            return PEPstrMODWrapper.submit_peptide_backup(sequence, options)

        # Extract PDB, remove preceding newline
        pdb = re.search(r"<pre>(?P<pdb>.+)</pre>", results_page.text, re.DOTALL).group("pdb")[1:]

        if cm().get("verbose"):
            logging.debug(f"Got result PDB from PEPstrMOD: {pdb}")

        return pdb

    @staticmethod
    def submit_peptide_backup(sequence: str, options: dict = None) -> str:
        """
        Submits a peptide sequence to the PEPstrMOD service hosted at a backup site (old one),
        optionally overriding the settings in the config.json, and returns the PDB as a str.
        This might take a few minutes.

        :param sequence: The peptide sequence in single letter format
        :param options: Any optional override settings
        :return: The PDB as a str
        """
        # TODO: Find out hardcoded values for all environment options
        form_data = {
            "MAX_FILE_SIZE": str(1000000),
            "seq": sequence,
            "email": cm().get("pepstrmod_config.email"),
            "time": cm().get("pepstrmod_config.simulation_time"),
            "env": cm().get("pepstrmod_config.environment"),
            "topol": "yes" if cm().get("pepstrmod_config.download_topology") else "no",
            "cluster": "yes" if cm().get("pepstrmod_config.cluster_analysis") else "no",
            "traj": "yes" if cm().get("pepstrmod_config.download_trajectory") else "no",
            "graph": "yes" if cm().get("pepstrmod_config.do_energy_rms_graph") else "no",
            "send": "Submit sequence for prediction",
        }

        if options:
            for k, v in options.items():
                form_data[k] = v

        # TODO: Do this by hand, so that the dependency can be dropped
        multipart_data = MultipartEncoder(
            fields={
                "MAX_FILE_SIZE": form_data["MAX_FILE_SIZE"],
                "seq": form_data["seq"],
                "email": form_data["email"],
                "time": form_data["time"],
                "env": form_data["env"],
                "topol": form_data["topol"],
                "cluster": form_data["cluster"],
                "traj": form_data["traj"],
                "graph": form_data["graph"],
                "send": form_data["send"],
            }
        )

        response = requests.post(url=PEPstrMODWrapper.ENDPOINT_REDUNDANCY, data=multipart_data,
                                 headers={'Content-Type': multipart_data.content_type})

        run_id = None
        # TODO: Think of better regex?
        pattern = re.compile(
            r'(https?://(www\.)?osddlinux\.osdd\.net/raghava/pepstrmod/results/index\.php\?ran=(?P<runid>\d+))')
        if matched_url := re.search(pattern, response.text):
            run_id = matched_url.group("runid")
            logging.debug(f"Identified PEPstrMOD run id {run_id}")
        else:
            raise ValueError("Could not identify where the PEPstrMOD results (run id) will be posted. The results "
                             "endpoint has likely changed, please create an issue for this on the GitHub repository.")

        pdb_url = f"http://osddlinux.osdd.net/tmp/pepstrmod/pep_{run_id}/pepstr{run_id}.pdb"
        logging.debug(f"Using PEPstrMOD endpoint {pdb_url}")

        # TODO: Find better marker for the results being ready?
        results_page = requests.get(pdb_url, timeout=300)
        interval = 60
        wait_limit = 30
        wait_count = 0
        while results_page.status_code == 404:
            if wait_count > wait_limit:
                # redundancy timed out, tertiary structure prediction currently not possible
                logging.error("Cannot predict tertiary peptide structure even when using backup service.")
                raise RuntimeError("Cannot predict tertiary peptide structure even when using backup service.")
            # Wait for results to be ready
            wait_count += 1
            time.sleep(interval)
            results_page = requests.get(pdb_url, allow_redirects=True, stream=False, timeout=300)

        # Extract PDB, remove preceding newline
        pdb = str(results_page.text)

        if cm().get("verbose"):
            logging.debug(f"Got result PDB from PEPstrMOD: {pdb}")

        return pdb

