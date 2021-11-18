###############
##   SETUP   ##
############### 
localrules: all

# Output file setup
output_files = ["early.csv", "late.csv", "seed.txt"]


# Parameters to vary
gens = [60]
Ks = [2750]
initPops = [2750]
iCorins = [0, 0.01, 0.05, 0.25] # expand
iEDNRBs = [0, 0.01, 0.05, 0.25] # expand
initOpts = [0.13] # Corresponds to NE Utah/ NW CO
finalOpts = [0.876] # Corresponds to NE Utah/ NW CO
selections = [0.6562] # Corresponding to a 5%, 7%, and 10% weekly survival penalty with 16 weeks of winter
H2s = [0.64]
offsprings = [4, 15] # 3.5 and 4 are close to what I've been working with, 15 is from James 1969
replicates = list(range(1,31)) # second number = 1 more than the desired number of replicates



# Filename pattern for the output files from each sim
res_pattern_additive_constant_K = "slim_results_final/additive_constantK/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_H2{H2}_lambda{offspring}_rep{replicate}_{output_file}"

# And expand all varying parameters to create the list of desired files
sim_results_additive_constantK = expand(res_pattern_additive_constant_K,                                
                                        gen = gens,
                                        K = Ks,
                                        initPop = initPops,
                                        iCorin = iCorins,
                                        iEDNRB = iEDNRBs,
                                        initOpt = initOpts,
                                        finalOpt = finalOpts,
                                        selection = selections,
                                        H2 = H2s,
                                        offspring = offsprings,
                                        replicate = replicates,
                                        output_file = output_files)


# Filename pattern for the output files from each sim
res_pattern_recessive_constant_K = "slim_results_final/recessive_constantK/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_H2{H2}_lambda{offspring}_rep{replicate}_{output_file}"

# And expand all varying parameters to create the list of desired files
sim_results_recessive_constantK = expand(res_pattern_recessive_constant_K,                                
                                        gen = gens,
                                        K = Ks,
                                        initPop = initPops,
                                        iCorin = iCorins,
                                        iEDNRB = iEDNRBs,
                                        initOpt = initOpts,
                                        finalOpt = finalOpts,
                                        selection = selections,
                                        H2 = H2s,
                                        offspring = offsprings,
                                        replicate = replicates,
                                        output_file = output_files)

# Filename pattern for output files from varying K sims
res_pattern_additive_varyK = "slim_results_final/additive_varyK/gens{gen}_minK{minK}_maxK{maxK}_period{period}_startK{startK}_initDec{initDec}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_H2{H2}_lambda{offspring}_rep{replicate}_{output_file}"

sim_results_additive_varyK = expand(res_pattern_additive_varyK,                                
                                    gen = gens,
                                    minK = [500],
                                    maxK = [5000],
                                    period = [9],
                                    startK = [2750],
                                    initDec = ["T","F"],
                                    iCorin = iCorins,
                                    iEDNRB = iEDNRBs,
                                    initOpt = [0.13],
                                    finalOpt = [0.876],
                                    selection = selections,
                                    H2 = H2s,
                                    offspring = offsprings,
                                    replicate = replicates,
                                    output_file = output_files)

###############
## Pipeline  ##
###############

# Pretty simple: one rule to run a single sim


rule all:
    input:
        sim_results_additive_constantK,
        sim_results_recessive_constantK,
        sim_results_additive_varyK

rule sim_additive_constantK:
    input:
    output:
        "slim_results_final/additive_constantK/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_H2{H2}_lambda{offspring}_rep{replicate}_early.csv",
        "slim_results_final/additive_constantK/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_H2{H2}_lambda{offspring}_rep{replicate}_late.csv",
        "slim_results_final/additive_constantK/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_H2{H2}_lambda{offspring}_rep{replicate}_seed.txt"
    log:
        "logs/slim/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_H2{H2}_lambda{offspring}_rep{replicate}.log"
    resources:
        cpus=1
    shell:
        """
        slim -d generations={wildcards.gen} -d K={wildcards.K} -d initPopSize={wildcards.initPop} -d initFreqCorin={wildcards.iCorin} -d initFreqEDNRB={wildcards.iEDNRB} -d initOptPheno={wildcards.initOpt} -d finalOptPheno={wildcards.finalOpt} -d fitFuncWidth={wildcards.selection} -d h2={wildcards.H2} -d offspringPoisLambda={wildcards.offspring}  -d replicate={wildcards.replicate} -d "outputDir='slim_results_final/additive_constantK'" src/slim_simulations/constant_K_additive_CLI.slim > {log} 
        """

rule sim_recessive_constantK:
    input:
    output:
        "slim_results_final/recessive_constantK/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_H2{H2}_lambda{offspring}_rep{replicate}_early.csv",
        "slim_results_final/recessive_constantK/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_H2{H2}_lambda{offspring}_rep{replicate}_late.csv",
        "slim_results_final/recessive_constantK/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_H2{H2}_lambda{offspring}_rep{replicate}_seed.txt"
    log:
        "logs/slim/recessive/gens{gen}_K{K}_initPop{initPop}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_H2{H2}_lambda{offspring}_rep{replicate}.log"
    resources:
        cpus=1
    shell:
        """
        slim -d generations={wildcards.gen} -d K={wildcards.K} -d initPopSize={wildcards.initPop} -d initFreqCorin={wildcards.iCorin} -d initFreqEDNRB={wildcards.iEDNRB} -d initOptPheno={wildcards.initOpt} -d finalOptPheno={wildcards.finalOpt} -d fitFuncWidth={wildcards.selection} -d h2={wildcards.H2} -d offspringPoisLambda={wildcards.offspring}  -d replicate={wildcards.replicate} -d "outputDir='slim_results_final/recessive_constantK'" src/slim_simulations/constant_K_recessive_CLI.slim > {log} 
        """

rule sim_additive_varyK:
    input:
    output:
        "slim_results_final/additive_varyK/gens{gen}_minK{minK}_maxK{maxK}_period{period}_startK{startK}_initDec{initDec}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_H2{H2}_lambda{offspring}_rep{replicate}_early.csv",
        "slim_results_final/additive_varyK/gens{gen}_minK{minK}_maxK{maxK}_period{period}_startK{startK}_initDec{initDec}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_H2{H2}_lambda{offspring}_rep{replicate}_late.csv",
        "slim_results_final/additive_varyK/gens{gen}_minK{minK}_maxK{maxK}_period{period}_startK{startK}_initDec{initDec}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_H2{H2}_lambda{offspring}_rep{replicate}_seed.txt"
    log:
        "logs/slim/additive_varyK/gens{gen}_minK{minK}_maxK{maxK}_period{period}_startK{startK}_initDec{initDec}_iCorin{iCorin}_iEDNRB{iEDNRB}_initOpt{initOpt}_finalOpt{finalOpt}_sel{selection}_H2{H2}_lambda{offspring}_rep{replicate}_seed.log"
    resources:
        cpus=1
    shell:
        """
        slim -d generations={wildcards.gen} -d min_K={wildcards.minK} -d max_K={wildcards.maxK} -d period={wildcards.period} -d start_K={wildcards.startK} -d initial_dec={wildcards.initDec} -d initFreqCorin={wildcards.iCorin} -d initFreqEDNRB={wildcards.iEDNRB} -d initOptPheno={wildcards.initOpt} -d finalOptPheno={wildcards.finalOpt} -d fitFuncWidth={wildcards.selection} -d h2={wildcards.H2} -d offspringPoisLambda={wildcards.offspring}  -d replicate={wildcards.replicate} -d "outputDir='slim_results_final/additive_varyK'" src/slim_simulations/varying_K_additive.slim > {log} 
        """