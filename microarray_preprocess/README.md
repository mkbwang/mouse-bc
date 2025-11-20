# Preprocess raw microarray files

* `annotationdb.py`: Utility functions for querying the NCBI database and ENSEMBL database for gene names.
* `microarray_ensembl_annotation.py` and `microarray_ensembl_annotation.sh` search the ENSEMBL databases for gene names corresponding to each accession.
* `microarray_refseq_annotation.py` and `microarray_refseq_annotation.sh` search the NCBI databases for gene names corresponding to each accession.
* `microarray_preprocess.R` does the RMA normalization of microarray data and match accession IDs to gene names.
