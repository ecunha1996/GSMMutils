import sys
import sympy as sp

def get_bounds(name, parameters):
    function_to_call = getattr(sys.modules[__name__], name)
    return function_to_call(parameters)


def light(parameters):
    Ke = 11.5 * parameters['X'] * parameters['chlorophyll']
    Ex0 = parameters['Eo'] / (parameters['Lr'] * Ke) * (1 - sp.exp(-parameters['Lr'] * Ke)) * parameters['light_conversion_factor']
    ro = parameters['ro1'] * parameters['chlorophyll'] + parameters['ro0']
    # Ex = ro / (parameters['X'] * parameters['Lr']) * Ex0 * parameters['light_conversion_factor']
    return Ex0


def nitrate(parameters):
    return sp.Max(sp.N(parameters["VNmax"] * parameters["nitrate"] / (parameters["KNm"] + parameters["nitrate"]) * (1 - parameters["q"])),0)


def starch_consumption(parameters):
    parameters["Ks"] = 0.034
    parameters["vmax"] = 0.066 / 48660.195 * 1000
    return sp.Max(0, sp.N(parameters["vmax"] * sp.Max(0, parameters["starch_concentration"]) / (sp.Max(0, parameters["starch_concentration"]) + parameters["Ks"])))


def starch_production(parameters):
    light_factor = parameters["Ex"] ** parameters["hill_coeff_starch"] / (parameters["Kstl"] ** parameters["hill_coeff_starch"] + parameters["Ex"] ** parameters["hill_coeff_starch"])
    return sp.Max(0, sp.N(parameters["maximum_starch_production"] * (1 - parameters["z"]) * light_factor))


def carotene(parameters):
    def phi(x, rs):
        return 1 / (1 + sp.exp(-rs * x))
    v_car_gen = parameters["v_car_max"] * (parameters["Ex"] ** parameters["l"]) / ((parameters["ExA"] ** parameters["l"]) + (parameters["Ex"] ** parameters["l"]))
    vcar = v_car_gen * phi(parameters["a1"] * parameters["Ex"] + parameters["a0"] - parameters["a2"] * parameters["nitrogen_mass_quota"] + parameters["a3"] * parameters["phosphate_mass_quota"] - parameters["a4"] * parameters["aeration"], parameters["smoothing_factor"])
    # vcar = v_car_gen * phi(parameters["a1"] * parameters["Ex"]- parameters["nitrogen_mass_quota"], parameters["smoothing_factor"]) * phi(parameters["a3"] * parameters["phosphate_mass_quota"], parameters["smoothing_factor"]) * phi(parameters["a4"] * parameters["aeration"], parameters["smoothing_factor"])
    # vcar = v_car_gen * phi(parameters["a1"] * parameters["Ex"] + parameters["a0"] - parameters["a2"] * parameters["nitrogen_mass_quota"] - parameters["a4"] * parameters["aeration"], parameters["smoothing_factor"]) * phi(parameters["a3"] * parameters["p_quota"], parameters["smoothing_factor_carotene"])
    # phosphate_factor = ((parameters["p_quota"] - parameters["wPmin"]) / (parameters["wPopt"] - parameters["wPmin"])) + parameters["c0"]
    # vcar = v_car_gen * phi(parameters["a1"] * parameters["Ex"] + parameters["a0"] - parameters["a2"] * parameters["nitrogen_mass_quota"], parameters["smoothing_factor"]) * phosphate_factor
    return sp.Max(0, sp.N(vcar))

def lutein(parameters):
    def phi(x, rs):
        return 1 / (1 + sp.exp(-rs * x))
    v_lut_gen = parameters["v_lut_max"] * (parameters["Ex"] ** parameters["l"]) / ((parameters["ExA"] ** parameters["l"]) + (parameters["Ex"] ** parameters["l"]))
    vlut = v_lut_gen * phi(parameters["a1_lut"] * parameters["Ex"] + parameters["a0_lut"] + parameters["a3_lut"] * parameters["phosphate_mass_quota"] - parameters["a4_lut"] * parameters["aeration"], parameters["smoothing_factor_lut"])
    return sp.Max(0, sp.N(vlut))


def chlorophyll(parameters):
    phosphate_factor = ((parameters["p_quota"] - parameters["wPmin"]) / (parameters["wPopt"] - parameters["wPmin"])) + parameters["c0"]
    def gamma(E):
        if E > parameters["Esat"]: E = parameters["Esat"]
        return parameters["ymax"] * (parameters["KEchl"] / (E + parameters["KEchl"]))

    def sum_chl(yE):
        return (yE - parameters["chlorophyll"] / parameters["nitrogen_mass_quota"]) * phosphate_factor
    return sum_chl(gamma(parameters["Ex"]))


def tag(parameters):
    return sp.Max(sp.N(parameters["maximum_tag_production"] / parameters["F"] / 904.78 * parameters["n"]), 0)


def glycerol(parameters):
    max_production = (1e-5 * parameters["nacl"] ** 2 + 0.002 * parameters["nacl"] + 0.112) / parameters["X"] * (1 - parameters["glycerol"] / parameters["wgly_max"])
    return sp.Max(0, sp.N(max_production))

def co2(parameters):
    return sp.Max(0, sp.N(parameters["vco2max"] * (1 - parameters["z"])))

