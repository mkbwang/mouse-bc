from annotationdb import ensembl_search
import argparse
import time
import numpy as np
import pandas as pd
import os


if __name__ == "__main__":

    folder="/nfs/turbo/sph-ligen/wangmk/mouse_bc/data/microarray"
    probe_metadata = pd.read_csv(os.path.join(folder, "probe_metadata.csv"))
    output_file = os.path.join(folder, "probe_ensembl_annotation.txt")

    n_probes = len(probe_metadata)

    db_names = probe_metadata["DB"].to_numpy()
    accessions = probe_metadata["accession"].to_numpy()

    for id in range(n_probes):
        dbname = db_names[id]
        accession = accessions[id]
        if dbname == "ENSEMBL":
            genename = ensembl_search(accession)
            time.sleep(1.5)
            output=f"{accession},{genename}\n"
            with open(output_file, 'a') as file:
                file.write(output)


