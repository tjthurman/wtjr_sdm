Dosage-env correlation README. 

This analysis does not require a conda environment, only a working version of Python (we used version 3.6.11) and the Numpy and Pandas libraries (we used versions 1.19.1 and 1.1.2, respectively). It also requires you to have run the genome-wide PCA code (in folder `02_genome-wide_PCA`), as some of the input files for this analysis are output by that pipeline. 

After running the PCA code, run the code in the `get_pcangsd_sites.sh` shell script to generate a file of the scaffolds and positions of the sites that went into the PCAngsd analysis. Then, run the `extract_color_gene_dosages.py` script to extract genotype dosages estimated by PCAngsd at the color-associated sites. 



These files show the code used in our analysis, but they will need to be modified (mostly by modifying paths) to run on a different computer. 
