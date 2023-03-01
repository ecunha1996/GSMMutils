from itertools import combinations_with_replacement, product

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from ExpAlgae.graphics.plot import lineplot
from ExpAlgae.model.COBRAmodel import MyModel
from os.path import join
from cobra.flux_analysis import pfba
from joblib import delayed, Parallel
from mewpy.simulation import get_simulator
import seaborn as sns

def convert_light_uptake(umol_m2_s):
    return umol_m2_s*5.7

def convert_light_uptake_rev(mmol_gDW_d):
    return mmol_gDW_d/5.7

def simulate_light_sources(model_copy: MyModel):
    light_sources = [reaction.id for reaction in model_copy.reactions if "PRISM" in reaction.id]
    results = {}
    for light_reaction in light_sources:
        with model_copy:
            model_copy.reactions.e_Biomass__cytop.bounds = (0.2, 0.2)
            model_copy.objective = "EX_C00205__dra"
            model_copy.reactions.get_by_id(light_reaction).bounds = (0, 1000)
            for other_reaction in light_sources:
                if other_reaction != light_reaction:
                    model_copy.reactions.get_by_id(other_reaction).bounds = [0, 0]
            sol = model_copy.optimize()
            results[light_reaction] = {'biomass': sol["e_Biomass__cytop"], "light_reaction": sol[light_reaction]}
    sns.barplot(list(key.split("__")[0] for key in results.keys()), list(val['light_reaction'] for val in results.values()))
    plt.show()


def run_simulation(model_copy, i):
    growth_rate =  [0, 0]
    with model_copy:
        constraint = model_copy.problem.Constraint(
            model_copy.reactions.R09542_hn437__chlo.flux_expression + model_copy.reactions.R09542_hn680__chlo.flux_expression,
            lb=i,
            ub=i)
        model_copy.add_cons_vars(constraint)
        # sol = model_copy.optimize()
        try:
            sol = pfba(model_copy)
        except:
            sol = model_copy.optimize()
        growth_rate[0] = sol["e_Biomass__cytop"]
        print(sol['R09503_hn438__lum'])
        print(sol['R09503_hn673__lum'])
    with model_copy:
        constraint = model_copy.problem.Constraint(
            model_copy.reactions.R09503_hn438__lum.flux_expression + model_copy.reactions.R09503_hn673__lum.flux_expression,
            lb=i,
            ub=i)
        model_copy.add_cons_vars(constraint)
        # sol = model_copy.optimize()
        try:
            sol = pfba(model_copy)
        except:
            sol = model_copy.optimize()
        growth_rate[1] = sol["e_Biomass__cytop"]
    return i, growth_rate


def main(model):
    light_range = np.arange(0, 4, 0.1)
    # for reaction in model.reactions:
    #     if reaction.lower_bound < -100:
    #         reaction.lower_bound = -10000
    #     if reaction.upper_bound > 100:
    #         reaction.upper_bound = 10000
    model.exchanges.EX_C00011__dra.bounds = (-8.21/24, 1000)
    model.reactions.ATPm__cytop.bounds = (2.85/24, 2.85/24)
    model.exchanges.EX_C00205__dra.bounds = (-1000, 1000)
    model.reactions.e_Pigment__chlo.bounds = (0, 10000)
    for reaction in model.reactions:
        if "trial" in reaction.id:
            reaction.bounds = (0, 0)
    result = Parallel(n_jobs=6)(delayed(run_simulation)(model.copy(), e) for e in light_range)
    res = [e[1] for e in result]
    light = [e[0] for e in result]
    growth_rate_ps2 = [e[1] for e in res]
    growth_rate_ps1 = [e[0] for e in res]
    # sns.lineplot([e for e in light_range], growth_rate_ps2).set(title='PSII')
    # plt.show()
    lineplot([e for e in light_range], growth_rate_ps2, title='PSII')
    lineplot([e for e in light_range], growth_rate_ps1, title='PSI')

    #create 3-surface plot
    # fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    # surf = ax.plot_surface(light_range, Y, Z, cmap=cm.coolwarm,
    #                        linewidth=0, antialiased=False)

    # experimental_x = [50, 200, 500, 1000, 1500]
    # experimental_y = [np.log(2)/198*24, np.log(2)/28*24, np.log(2)/7.5*24, np.log(2)/7.2*24, np.log(2)/6.5*24]
    # sns.scatterplot(experimental_x, experimental_y)


def run_simulation_2(model_copy, fluxes):
    growth_rate =  0
    with model_copy:
        constraint = model_copy.problem.Constraint(
            model_copy.reactions.R09542_hn437__chlo.flux_expression + model_copy.reactions.R09542_hn680__chlo.flux_expression,
            lb=fluxes[0],
            ub=fluxes[0])
        model_copy.add_cons_vars(constraint)
        constraint = model_copy.problem.Constraint(
            model_copy.reactions.R09503_hn438__lum.flux_expression + model_copy.reactions.R09503_hn673__lum.flux_expression,
            lb=fluxes[1],
            ub=fluxes[1])
        model_copy.add_cons_vars(constraint)
        # sol = model_copy.optimize()
        try:
            sol = pfba(model_copy)
        except:
            sol = model_copy.optimize()
        growth_rate = sol["e_Biomass__cytop"]
    return fluxes, growth_rate


def main_ps(data_directory, model):
    light_range_ps1 = np.arange(0, 12, 0.1)
    light_range_ps2 = np.arange(0, 20, 0.1)
    model.exchanges.EX_C00011__dra.bounds = (-8.21 / 24, 1000)
    model.reactions.ATPm__cytop.bounds = (2.85 / 24, 2.85 / 24)
    model.exchanges.EX_C00205__dra.bounds = (-1000, 1000)
    model.reactions.e_Pigment__chlo.bounds = (0, 1000)
    for reaction in model.reactions:
        if "trial" in reaction.id:
            reaction.bounds = (0, 0)
    result = Parallel(n_jobs=7)(delayed(run_simulation_2)(model.copy(), e) for e in list(product(light_range_ps2.tolist(), light_range_ps1.tolist())))
    # create 3-surface plot
    df = pd.DataFrame(data=list(zip([e[0][0] for e in result], [e[0][1] for e in result], [e[1] for e in result])), columns=['PSI', 'PSII', 'growth_rate'])
    df.to_csv(join(data_directory, "ps.csv"), index=False)
    df = pd.read_csv(join(data_directory, "ps.csv"))
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf=ax.plot_trisurf( df['PSII'], df['PSI'], df['growth_rate'], cmap=plt.cm.jet, linewidth=0.2)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.savefig(join(data_directory, "ps.png"))


if __name__ == '__main__':
    data_directory = "../data"
    model = MyModel(join(data_directory, "models/model_with_trials.xml"), "e_Biomass__cytop")
    # main(model)
    main_ps(data_directory, model)
