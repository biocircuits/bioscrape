import bioscrape as bs
import bioscrape.lineage as bs_lineage
import math
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

def make_growth_model(rxns, species, params, initial_conditions):
    '''
    Adds some scaffolding to a model to put it in an exponentially-growing cell 
    with growth rate params['k_gamma']. Parameters are of the same type and 
    specification as used by bioscrape.types.Model.
    '''
    m = bs_lineage.LineageModel(reactions = rxns, 
                              parameters = params,
                              species = species,
                              initial_condition_dict = initial_conditions)

    # Multiplicative growth with rate set by k_gamma, but don't grow if there's 
    # no DNA present. 
    m.create_volume_rule("ode", {"equation": "volume * _k_gamma * Heaviside(DNA-1)"})

    # Divide when volume doubles (to 2), splitting all species binomially, with a small amount
    # of noise in volume partitioning.
    vsplit = bs_lineage.LineageVolumeSplitter(m, 
              options = {"default": "binomial"},
              partition_noise = 0.05)
    division_vol = 2
    m.create_division_rule("volume", {"threshold":division_vol}, vsplit)

    m.py_initialize()
    return m

def make_trivial_model(params, initial_conditions, lineage = False):
    rxns = [
        (['DNA'], ['DNA', 'DNA'], 'massaction', {'k': params['k_alpha']})
    ]
    if not lineage:
        rxns.append((['DNA'], [], 'massaction', {'k': params['k_gamma']}))
        return bs.types.Model(species = ["DNA"], parameters = params, reactions = rxns, 
                initial_condition_dict = initial_conditions)
    else:
        m = make_growth_model(rxns, ["DNA"], params, initial_conditions)
        return m

def make_dummy_model(params, initial_conditions, lineage = False):
    rxns = [
        ([], ['R'], 'massaction', {'k': params['k_alpha']}),
        (['R', 'DNA'], ['DNA', 'DNA'], 'massaction', {'k':params['k_fast']})
    ]
    if not lineage:
        rxns.append((['DNA'], [], 'massaction', {'k':params['k_gamma']}))
        return bs.types.Model(species = ['DNA', 'R'], parameters = params, reactions = rxns,
                              initial_condition_dict = initial_conditions)
    else:
        m = make_growth_model(rxns, ["R", "DNA"], params, initial_conditions)
        return m

# These are the parameters reported in the original B&P paper, with division
# time set to 30 minutes. 
bp_default_params = {
        "k_II": 0.25,
        "k_L": 12.0,
        "k_mL": 4.3,
        "k_p": 4.3,
        "k_D": 5,
        "k_1": 0.15, # 1/nM
        "k_m1": 48,
        "k_2": 44,
        "k_m2": 0.085,
        "k_mC": 17,
        "k_I": 6,
        "k_gamma_I": .35,
        "k_gamma": math.log(2)/30 # 30 min division time
    }

# The same parameters, but nondimensionalized so that division time = 1
bp_nondim_params = {param: val/(bp_default_params['k_gamma']) for (param,val) in bp_default_params.items()}

def make_BP_model(params, initial_conditions, lineage = False):
    plasmid_species = ["DNA", "DNA_RIIs", "DNA_RIIL", "DNA_p", "DNA_II_Iu", "DNA_II_Is"]
    species = plasmid_species + ["RI"]

    reactions = [
        (["DNA"], ["DNA_RIIs"], "massaction", {"k": "k_II"}),
        (["DNA"], ["DNA", "RI"], "massaction", {"k": "k_I"}),
        (["DNA_RIIs"], ["DNA_RIIL"], "massaction", {"k": "k_L"}),
        (["DNA_RIIL"], ["DNA"], "massaction", {"k": "k_mL"}),
        (["DNA_RIIL"], ["DNA_p"], "massaction", {"k": "k_p"}),
        (["DNA_p"], ["DNA", "DNA"], "massaction", {"k": "k_D"}),
        (["DNA_RIIs", "RI"], ["DNA_II_Iu"], "massaction", {"k": "k_1"}),
        (["DNA_II_Iu"], ["DNA_RIIs", "RI"], "massaction", {"k": "k_m1"}),
        (["DNA_II_Iu"], ["DNA_II_Is"], "massaction", {"k": "k_2"}),
        (["DNA_II_Is"], ["DNA_II_Iu"], "massaction", {"k": "k_m2"}),
        (["DNA_II_Is"], ["DNA"], "massaction", {"k": "k_mC"}),
        (["RI"], [], "massaction", {"k": "k_gamma_I"})
    ]
    
    if not lineage:
        for s in species:
            reactions.append(([s], [], "massaction", {"k": "k_gamma"}))
        return bs.types.Model(species = species, parameters = params, reactions = reactions,
                              initial_condition_dict = initial_conditions)
    else:
        m = make_growth_model(reactions, species, params, initial_conditions)
        return m

def make_simplified_BP_model(params, initial_conditions, lineage = False):
    species = ["DNA", "DNAp", "R"]

    reactions = [
        (["DNA"], ["DNAp"], "massaction", {"k": "k_p"}),
        (["DNAp"], ["DNA", "DNA"], "massaction", {"k": "k_rep"}),
        (["DNA"], ["DNA", "R"], "massaction", {"k": "k_tx"}),
        (["DNAp", "R"], ["DNA"], "massaction", {"k": "k_I"}),
        (["R"], [], "massaction", {"k": "k_gamma_I"})
    ]
    
    if not lineage:
        for s in species:
            reactions.append(([s], [], "massaction", {"k": "k_gamma"}))
        return bs.types.Model(species = species, parameters = params, reactions = reactions,
                              initial_condition_dict = initial_conditions)
    else:
        m = make_growth_model(reactions, species, params, initial_conditions)
        return m



def main():
    default_trivial_params = {
        "k_gamma": 1,
        "k_alpha": 1
    }   

    default_dummy_params = {
        "k_gamma": 1,
        "k_alpha": 50,
        "k_fast": 100
    }

    # These are the parameters reported in the original B&P paper, with division
    # time set to 30 minutes. 
    bp_default_params = {
            "k_II": 0.25,
            "k_L": 12.0,
            "k_mL": 4.3,
            "k_p": 4.3,
            "k_D": 5,
            "k_1": 0.15, # 1/nM
            "k_m1": 48,
            "k_2": 44,
            "k_m2": 0.085,
            "k_mC": 17,
            "k_I": 6,
            "k_gamma_I": .35,
            "k_gamma": math.log(2)/30 # 30 min division time
        }

    # The same parameters, but nondimensionalized so that division time = 1
    bp_nondim_params = {param: val/(bp_default_params['k_gamma']) for (param,val) in bp_default_params.items()}

    # Simplified B&P model parameters, fit against data generated from BP model
    # (not shown)
    fit_simple_bp_params = {
        "k_rep": 14.156455,
        "k_tx": 252.232848,
        "k_I": 0.665226,
        "k_p": 36.603267,
        "k_gamma_I": 12.534647,
        "k_gamma": 1
    }

    all_distributions = dict()
    all_models = dict()

    default_init = {"DNA": 50}
    all_models['trivial']   = make_trivial_model(default_trivial_params, default_init, lineage = True)
    all_models['dummy']     = make_dummy_model(default_dummy_params, default_init, lineage = True)
    all_models['bp']        = make_BP_model(bp_nondim_params, default_init, lineage = True)
    all_models['simple_bp'] = make_simplified_BP_model(fit_simple_bp_params, default_init, lineage = True)

    timings = {name:[] for name, m in all_models.items()}
    Ns = [int(n) for n in np.geomspace(2, 2048, 10)]

    for model_name, model in all_models.items():
        if model_name == "trivial":
            continue
        for population_cap in Ns:
            seed = 42334
            np.random.seed(42334)
            bs.random.py_seed_random(42334)
            print("Simulating distribution of model " + model_name + f" with N = {population_cap}...")
            # The trivial model spreads more over time, so worth simulating longer
            # to get more distributional info.
            n_generations = 30 if model_name == "trivial" else 50
            temp_ts = np.linspace(0, n_generations, 100)#n_generations*100)
            sample_times = np.linspace(0, n_generations, 40) #<- Lineages doesn't like having too few
                                                            # sample times, don't use just one or two.

        #     all_distributions[model_name] = bs_lineage.py_SimulateCellLineage(temp_ts, [], initial_cell_count = 1, Model = model)
        #     [bs.simulator.py_simulate_model(temp_ts, model, stochastic = True) for i in range(2)]
            loop_start_time = datetime.now()
            _ = bs_lineage.py_SimulateTurbidostat(initial_cell_states = 1, 
                                                  timepoints = temp_ts,
                                                  sample_times = sample_times,
                                                  population_cap = population_cap,
                                                  Model = model, debug = False)
            loop_end_time = datetime.now()
            timings[model_name].append(loop_end_time - loop_start_time)
    print("Done simulating.")
    print("Plotting...")

    for model_name, model in all_models.items():
        if model_name == "trivial":
            continue
        seconds_per_cell = [timings[model_name][i].total_seconds()/Ns[i] for i in range(len(Ns))]
        plt.scatter(Ns, seconds_per_cell, label = model_name)
    plt.legend()
    plt.xlabel("Population Cap")
    plt.ylabel("Time/Cell (real-time)")
    plt.title("py_SimualateTurbidostat Performance")
    plt.show()

    print("Done.")

if __name__ == "__main__":
    main()


