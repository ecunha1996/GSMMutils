import os
from os.path import join
import sys
from dfba import DfbaModel, ExchangeFlux, KineticVariable
from plotly import graph_objects as go
from dfba.plot.plotly import *
sys.path.insert(0, "/home/src")

from ExpAlgae.experimental.ExpMatrix import *
from ExpAlgae.model.COBRAmodel import MyModel


os.chdir("/home/")

def read_model(data_directory):
    model = MyModel(join(data_directory, "model_with_biomass_trials.xml"), "e_Biomass__cytop")
    model.add_medium(join(data_directory, "media.xlsx"), "base_medium")
    model.reactions.e_Biomass_ht__cytop.bounds = (0, 0)
    model.reactions.R00019__chlo.bounds = (0, 0)
    model.reactions.R00019__mito.bounds = (0, 0)
    model.set_prism_reaction("PRISM_white_LED__extr")
    model.reactions.ATPm__cytop.bounds = (2.85, 2.85)
    return model


def create_dfba_model(fba_model, matrix, index):
    dfba_model = DfbaModel(fba_model)
    # dfba_model.solver_data.set_algorithm("direct")
    dfba_model.solver_data.set_display("none")
    dfba_model.solver_data.set_ode_method("ADAMS")
    X = KineticVariable("Biomass")
    nitrate = KineticVariable("Nitrate")
    phosphorus = KineticVariable("Phosphate")
    starch = KineticVariable("Starch")
    tag = KineticVariable("TAG")
    carotene = KineticVariable("Carotene")
    M_p = fba_model.metabolites.C00009__extr.formula_weight
    M_n = fba_model.metabolites.C00244__extr.formula_weight
    dfba_model.add_kinetic_variables([X, nitrate, phosphorus, starch, tag, carotene])

    mu = ExchangeFlux(f"e_Biomass_trial{index}__cytop")
    v_N = ExchangeFlux("EX_C00244__dra")
    v_P = ExchangeFlux("EX_C00009__dra")
    v_S = ExchangeFlux("DM_C00369__chlo")
    v_T = ExchangeFlux("DM_C00422__lip")
    v_C = ExchangeFlux("DM_C02094__chlo")

    dfba_model.add_exchange_fluxes([mu, v_N, v_P, v_S, v_T, v_C])

    dfba_model.add_rhs_expression("Biomass", mu * X)
    dfba_model.add_rhs_expression("Nitrate", v_N * X)
    dfba_model.add_rhs_expression("Phosphate", v_P * X)
    dfba_model.add_rhs_expression("Starch", v_S * X)
    dfba_model.add_rhs_expression("TAG", v_T * X)
    dfba_model.add_rhs_expression("Carotene", v_C * X)

    dfba_model.add_exchange_flux_lb("EX_C00009__dra", 0.39 * phosphorus / (0.0038 + phosphorus), phosphorus)

    dfba_model.add_initial_conditions(
        {
            "Biomass": matrix.matrix[index]["DW"][0],
            "Nitrate": matrix.conditions["[N] mmol"].loc[index],
            "Phosphate": matrix.conditions["[P] mmol"].loc[index],
        }
    )
    concentrations, trajectories = dfba_model.simulate(
        0.0, max(matrix.matrix[index].index.astype(int).tolist()), 1, [f"e_Biomass_trial{index}__cytop", "EX_C00009__dra", "DM_C02094__chlo", "DM_C00369__chlo", "DM_C00116__cytop", "DM_C00422__lip"]
    )
    experimental = matrix.matrix[index].index, matrix.matrix[index]["DW"].tolist()

    fig = plot_concentrations(concentrations, metabolites=["Phosphate"], right_y_axis_title=r"$\textrm{Phosphate} \left[ \textrm{mmol} \, "
                                                                                            r"\textrm{L}^{-1} \right]$")
    fig.add_trace(
        go.Scatter(
            x=experimental[0],
            y=experimental[1],
            mode="markers",
            name="Biomass (exp.)",
        ),
    )
    fig.update_layout(
        autosize=False,
        width=700,
        height=500)

    fig.write_image(f"/home/data/concentrations_{index}.png")

    fig = plot_trajectories(trajectories, x_axis_title=r"$\textrm{Time} \left[ \textrm{h} \right]$", y_axis_title=r"$\textrm{Flux} \left[ \textrm{mmol} \, "
                                                                                                                  r"\textrm{g}_{\textrm{DW}}^{-1} \, \textrm{h}^{-1} "
                                                                                                                  r"\right]$")
    fig.write_image(f"/home/data/trajectories_{index}.png")

    concentrations.to_csv(f"/home/data/concentrations_{index}.csv", index=True)
    trajectories.to_csv(f"/home/data/trajectories_{index}.csv", index=True)


def simulation_for_conditions(matrix):
    as_dict = matrix.conditions[["C00011"]].to_dict(orient='index')
    for index, condition in as_dict.items():
        fba_model = read_model("/home/data")
        fba_model.exchanges.EX_C00011__dra.bounds = (-1000, 1000)
        for reaction in fba_model.reactions:
            if "Biomass" in reaction.id and "EX_" not in reaction.id and reaction.id != f"e_Biomass_trial{index}__cytop":
                reaction.bounds = (0, 0)
        fba_model.reactions.get_by_id(f"e_Biomass_trial{index}__cytop").bounds = (0, 1000)
        fba_model.objective = f"e_Biomass_trial{index}__cytop"
        for met, lb in condition.items():
            lb = -lb if lb < 0 else lb
            fba_model.reactions.get_by_id("EX_" + met + "__dra").bounds = (round(-lb, 4), 1000)
        create_dfba_model(fba_model, matrix, index)



matrix = ExpMatrix("/home/data/Matriz- DCCR Dunaliella salina_new.xlsx")
matrix.conditions = "Resume"
simulation_for_conditions(matrix)



# experimental = [0, 2, 4, 6, 8, 10, 12], [0.126, 0.221, 0.371, 0.458, 0.532, 0.576, 0.581]

#
# experimental = [0, 48, 96, 144, 192, 240, 288, 336, 384], \
#                [0.137, 0.214, 0.373, 0.616,
# 0.895,
# 1.035,
# 1.064,
# 1.161,
# 1.173
# ]

# experimental = [0, 48, 96, 144, 192, 240, 288, 336, 384], \
#                [0.107,
# 0.182,
# 0.368,
# 0.482,
# 0.615,
# 0.704,
# 0.762,
# 0.824,
# 0.792
# ]

# experimental = [0, 2, 4, 6, 8, 10, 12, 14, 16], [0.137, 0.214, 0.373, 0.616,
# 0.895,
# 1.035,
# 1.064,
# 1.161,
# 1.173
# ]

# fig = plot_concentrations(concentrations, y= "Phosphate", experimental=experimental)



# concentrations = pd.read_csv("concentrations.csv")
# trajectories = pd.read_csv("trajectories.csv")   # , "Starch", "Carotene", "TAG"


"""Instructions:
1. run in the container terminal: python3 run_dfba.py
2. install pip install -U kaleido, plotly
3. change solver to cplex
"""


"""
Calculations:


Vmax (P):

Vmax = 220 umol/molC/min
Vmax = 220 * 6.5458315 * 60 /1000 = 86.4049758 (valor absurdo)
OR
Vmax = 7 amol/um^3/min
Vmax = 

OR
Vmax = 1fm/cell/min
Vmax = 1*10^-15 / 153*10^-12 = 0.39 mmol/g/h

Vmax = 28 pmol/cell/min
Vmax = 28*10^-12 / 153*10^-12 *60  = 0.011 mmol/g/h (control)
Vmax = 730*10^-12 / 153*10^-12 *60  =  0.28627451 mmol/g/h (P depleted)

Km (P):

Km = 3.8 umol/L = 0.0038 mmol/L

Km = 1.4 mol/L = 0.0014 mmol/L

"""