# chang to folder with all the scaffold GL files
/home/tt164677e/tim_beegfs/wtjr_genomics/results/angsd/beagle_GL_scaffold/all_wtjr_samples/by_scaffold

# Corin scaffold 342
zcat 342_GL.beagle.gz | grep -E '342_47124004|342_46991691|342_46991393|342_46989842|342_46983691|342_46972767|342_46970080|342_46966737|342_46958289|342_46957711|342_46879687|342_46871251|342_46870994' > ../corin_site_GLs.beagle

# Agouti scaffold 245
zcat 245_GL.beagle.gz | grep -E '245_24236852|245_24229551|245_24226798|245_24225226|245_24222426|245_24221492|245_24220190|245_24211740|245_24204540|245_24200684|245_24197147' > ../asip_site_GLs.beagle

# EDNRB scaffold 311
zcat 311_GL.beagle.gz | grep -E '311_3467105|311_3456263|311_3454568|311_3419309|311_3419220|311_3413288|311_3402184|311_3399522' > ../ednrb_site_GLs.beagle

# Scaffold 380 unknown gene
zcat 380_GL.beagle.gz | grep -E '380_35303819|380_35303804|380_35321971|380_35322089|380_35311674|380_35320825|380_35253666|380_35322019|380_35313330|380_35306181' > ../scaff380_site_GLs.beagle