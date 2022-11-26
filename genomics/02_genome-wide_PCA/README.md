The `angsd.yml` file describes a conda environment needed to perform the genome-wide PCA analysis. PCAngsd is not available from `conda`, it must be installed separately from http://www.popgen.dk/software/index.php/PCAngsd. We used PCAngsd v1.03. The snakemake file presents the code used in our analysis, but would need to be modified to run properly on a different computer (updating the paths to the bamlist file, reference genome, and PCAngsd).