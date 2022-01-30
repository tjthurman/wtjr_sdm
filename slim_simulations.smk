###############
##   SETUP   ##
############### 
localrules: all, compile_slim_early, compile_slim_late, compile_slim_seed

# Output file setup
output_files = ["early.csv", "late.csv", "seed.txt"]


# Parameters to vary
gens = [80]
Ks = [2750]
initPops = [2750]
iCorins = [0, 0.01, 0.05, 0.25] # expand
iEDNRBs = [0, 0.01, 0.05, 0.25] # expand
initOpts = [0.13] # Corresponds to NE Utah/ NW CO
finalOpts = [0.876] # Corresponds to NE Utah/ NW CO
selections = [0.7805, 0.6562, 0.5446] # Corresponding to a 5%, 7%, and 10% weekly survival penalty with 16 weeks of winter
VEs = [3.59206] # From previous SliM simulations
offsprings = [4, 15] # 15 is from James 1969
replicates = list(range(1,31)) # second number = 1 more than the desired number of replicates



# Filename pattern for the output files from each individual simulation 
# with the additive, constantK, two locus model
res_pattern_additive_consK_2locus = "results/slim_ind_sims/additive_constantK_2locus/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_Ve{VE}_lambda{offspring}_rep{replicate}_{output_file}"

# And expand all varying parameters to create the list of desired files
sim_results_additive_consK_2locus = expand(res_pattern_additive_consK_2locus,                                
                                        gen = gens,
                                        K = Ks,
                                        initPop = initPops,
                                        iCorin = iCorins,
                                        iEDNRB = iEDNRBs,
                                        initOpt = initOpts,
                                        finalOpt = finalOpts,
                                        selection = selections,
                                        VE = VEs,
                                        offspring = offsprings,
                                        replicate = replicates,
                                        output_file = output_files)


# Filename pattern for the output files from each individual simulation 
# with the recessive, constantK, two locus model
res_pattern_recessive_consK_2locus = "results/slim_ind_sims/recessive_constantK_2locus/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_Ve{VE}_lambda{offspring}_rep{replicate}_{output_file}"

# And expand all varying parameters to create the list of desired files
sim_results_recessive_consK_2locus = expand(res_pattern_recessive_consK_2locus,                                
                                        gen = gens,
                                        K = Ks,
                                        initPop = initPops,
                                        iCorin = iCorins,
                                        iEDNRB = iEDNRBs,
                                        initOpt = initOpts,
                                        finalOpt = finalOpts,
                                        selection = selections,
                                        VE = VEs,
                                        offspring = offsprings,
                                        replicate = replicates,
                                        output_file = output_files)

# Filename pattern for the output files from each individual simulation 
# with the additive, varying K, two locus model
res_pattern_additive_varyK_2locus = "results/slim_ind_sims/additive_varyK_2locus/gens{gen}_minK{minK}_maxK{maxK}_period{period}_startK{startK}_initDec{initDec}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_Ve{VE}_lambda{offspring}_rep{replicate}_{output_file}"

sim_results_additive_varyK_2locus = expand(res_pattern_additive_varyK_2locus,                                
                                    gen = gens,
                                    minK = [500],
                                    maxK = [5000],
                                    period = [9],
                                    startK = [2750],
                                    initDec = ["T","F"],
                                    iCorin = iCorins,
                                    iEDNRB = iEDNRBs,
                                    initOpt = initOpts,
                                    finalOpt = finalOpts,
                                    selection = selections,
                                    VE = VEs,
                                    offspring = offsprings,
                                    replicate = replicates,
                                    output_file = output_files)
                                

# Filename pattern for the output files from each individual simulation 
# with the recessive, constant K, one locus model
res_pattern_recessive_consK_1locus = "results/slim_ind_sims/recessive_constantK_1locus/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_Ve{VE}_lambda{offspring}_rep{replicate}_{output_file}"

# And expand all varying parameters to create the list of desired files
sim_results_recessive_consK_1locus = expand(res_pattern_recessive_consK_1locus,                                
                                        gen = gens,
                                        K = Ks,
                                        initPop = initPops,
                                        iCorin = iCorins,
                                        initOpt = initOpts,
                                        finalOpt = finalOpts,
                                        selection = selections,
                                        VE = VEs,
                                        offspring = offsprings,
                                        replicate = replicates,
                                        output_file = output_files)



###############
## Pipeline  ##
###############

rule all:
    input: # slim simulation summaries
        # Early files:
        "results/slim_summaries/additive_constantK_2locus_early.csv",
        "results/slim_summaries/additive_varyK_2locus_early.csv",
        "results/slim_summaries/recessive_constantK_2locus_early.csv",
        "results/slim_summaries/recessive_constantK_1locus_early.csv",
        # Late files
        "results/slim_summaries/additive_constantK_2locus_late.csv",
        "results/slim_summaries/additive_varyK_2locus_late.csv",
        "results/slim_summaries/recessive_constantK_2locus_late.csv",
        "results/slim_summaries/recessive_constantK_1locus_late.csv",
        # seed files
        "results/slim_summaries/additive_constantK_2locus_seed.csv",
        "results/slim_summaries/additive_varyK_2locus_seed.csv",
        "results/slim_summaries/recessive_constantK_2locus_seed.csv",
        "results/slim_summaries/recessive_constantK_1locus_seed.csv"

rule sim_additive_constantK_2locus:
    input:
    output:
        "results/slim_ind_sims/additive_constantK_2locus/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_Ve{VE}_lambda{offspring}_rep{replicate}_early.csv",
        "results/slim_ind_sims/additive_constantK_2locus/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_Ve{VE}_lambda{offspring}_rep{replicate}_late.csv",
        "results/slim_ind_sims/additive_constantK_2locus/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_Ve{VE}_lambda{offspring}_rep{replicate}_seed.txt"
    log:
        "logs/slim/additive_consK_2locus/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_Ve{VE}_lambda{offspring}_rep{replicate}.log"
    resources:
        cpus=1
    shell:
        """
        slim -d generations={wildcards.gen} -d K={wildcards.K} -d initPopSize={wildcards.initPop} -d initFreqCorin={wildcards.iCorin} -d initFreqEDNRB={wildcards.iEDNRB} -d initOptPheno={wildcards.initOpt} -d finalOptPheno={wildcards.finalOpt} -d fitFuncWidth={wildcards.selection} -d V_E={wildcards.VE} -d offspringPoisLambda={wildcards.offspring}  -d replicate={wildcards.replicate} -d "outputDir='results/slim_ind_sims/additive_constantK_2locus'" src/slim_simulations/constant_K_additive_2locus_consVE_CLI.slim 2> {log} 
        """

rule sim_recessive_constantK_2locus:
    input:
    output:
        "results/slim_ind_sims/recessive_constantK_2locus/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_Ve{VE}_lambda{offspring}_rep{replicate}_early.csv",
        "results/slim_ind_sims/recessive_constantK_2locus/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_Ve{VE}_lambda{offspring}_rep{replicate}_late.csv",
        "results/slim_ind_sims/recessive_constantK_2locus/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_Ve{VE}_lambda{offspring}_rep{replicate}_seed.txt"
    log:
        "logs/slim/recessive_consK_2locus/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_Ve{VE}_lambda{offspring}_rep{replicate}.log"
    resources:
        cpus=1
    shell:
        """
        slim -d generations={wildcards.gen} -d K={wildcards.K} -d initPopSize={wildcards.initPop} -d initFreqCorin={wildcards.iCorin} -d initFreqEDNRB={wildcards.iEDNRB} -d initOptPheno={wildcards.initOpt} -d finalOptPheno={wildcards.finalOpt} -d fitFuncWidth={wildcards.selection} -d V_E={wildcards.VE} -d offspringPoisLambda={wildcards.offspring}  -d replicate={wildcards.replicate} -d "outputDir='results/slim_ind_sims/recessive_constantK_2locus'" src/slim_simulations/constant_K_recessive_2locus_consVE_CLI.slim 2> {log} 
        """

rule sim_additive_varyK_2locus:
    input:
    output:
        "results/slim_ind_sims/additive_varyK_2locus/gens{gen}_minK{minK}_maxK{maxK}_period{period}_startK{startK}_initDec{initDec}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_Ve{VE}_lambda{offspring}_rep{replicate}_early.csv",
        "results/slim_ind_sims/additive_varyK_2locus/gens{gen}_minK{minK}_maxK{maxK}_period{period}_startK{startK}_initDec{initDec}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_Ve{VE}_lambda{offspring}_rep{replicate}_late.csv",
        "results/slim_ind_sims/additive_varyK_2locus/gens{gen}_minK{minK}_maxK{maxK}_period{period}_startK{startK}_initDec{initDec}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_Ve{VE}_lambda{offspring}_rep{replicate}_seed.txt"
    log:
        "logs/slim/additive_varyK_2locus/gens{gen}_minK{minK}_maxK{maxK}_period{period}_startK{startK}_initDec{initDec}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_Ve{VE}_lambda{offspring}_rep{replicate}_seed.log"
    resources:
        cpus=1
    shell:
        """
        slim -d generations={wildcards.gen} -d min_K={wildcards.minK} -d max_K={wildcards.maxK} -d period={wildcards.period} -d start_K={wildcards.startK} -d initial_dec={wildcards.initDec} -d initFreqCorin={wildcards.iCorin} -d initFreqEDNRB={wildcards.iEDNRB} -d initOptPheno={wildcards.initOpt} -d finalOptPheno={wildcards.finalOpt} -d fitFuncWidth={wildcards.selection} -d V_E={wildcards.VE} -d offspringPoisLambda={wildcards.offspring}  -d replicate={wildcards.replicate} -d "outputDir='results/slim_ind_sims/additive_varyK_2locus'" src/slim_simulations/varying_K_additive_2locus_consVE_CLI.slim 2> {log} 
        """
    

rule sim_recessive_constantK_1locus:
    input:
    output:
        "results/slim_ind_sims/recessive_constantK_1locus/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_Ve{VE}_lambda{offspring}_rep{replicate}_early.csv",
        "results/slim_ind_sims/recessive_constantK_1locus/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_Ve{VE}_lambda{offspring}_rep{replicate}_late.csv",
        "results/slim_ind_sims/recessive_constantK_1locus/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_Ve{VE}_lambda{offspring}_rep{replicate}_seed.txt"
    log:
        "logs/slim/recessive_consK_1locus/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_Ve{VE}_lambda{offspring}_rep{replicate}.log"
    resources:
        cpus=1
    shell:
        """
        slim -d generations={wildcards.gen} -d K={wildcards.K} -d initPopSize={wildcards.initPop} -d initFreqCorin={wildcards.iCorin} -d initOptPheno={wildcards.initOpt} -d finalOptPheno={wildcards.finalOpt} -d fitFuncWidth={wildcards.selection} -d V_E={wildcards.VE} -d offspringPoisLambda={wildcards.offspring}  -d replicate={wildcards.replicate} -d "outputDir='results/slim_ind_sims/recessive_constantK_1locus'" src/slim_simulations/constant_K_recessive_1locus_consVE_CLI.slim 2> {log} 
        """

# Compile the many inidividual text files into single .csvs for later use. 
rule compile_slim_early:
    input:
        sim_results_additive_consK_2locus,
        sim_results_recessive_consK_2locus,
        sim_results_additive_varyK_2locus,
        sim_results_recessive_consK_1locus
    output:
        "results/slim_summaries/additive_constantK_2locus_early.csv",
        "results/slim_summaries/additive_varyK_2locus_early.csv",
        "results/slim_summaries/recessive_constantK_2locus_early.csv",
        "results/slim_summaries/recessive_constantK_1locus_early.csv"
    resources:
        cpus=1
    shell:
        """
        Rscript src/compile_slim_early.R
        """

rule compile_slim_late:
    input:
        sim_results_additive_consK_2locus,
        sim_results_recessive_consK_2locus,
        sim_results_additive_varyK_2locus,
        sim_results_recessive_consK_1locus
    output:
        "results/slim_summaries/additive_constantK_2locus_late.csv",
        "results/slim_summaries/additive_varyK_2locus_late.csv",
        "results/slim_summaries/recessive_constantK_2locus_late.csv",
        "results/slim_summaries/recessive_constantK_1locus_late.csv"
    resources:
        cpus=1
    shell:
        """
        Rscript src/compile_slim_late.R
        """

rule compile_slim_seed:
    input:
        sim_results_additive_consK_2locus,
        sim_results_recessive_consK_2locus,
        sim_results_additive_varyK_2locus,
        sim_results_recessive_consK_1locus
    output:
        "results/slim_summaries/additive_constantK_2locus_seed.csv",
        "results/slim_summaries/additive_varyK_2locus_seed.csv",
        "results/slim_summaries/recessive_constantK_2locus_seed.csv",
        "results/slim_summaries/recessive_constantK_1locus_seed.csv"
    resources:
        cpus=1
    shell:
        """
        Rscript src/compile_slim_seed.R
        """