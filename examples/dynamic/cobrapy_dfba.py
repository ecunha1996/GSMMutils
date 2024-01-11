from os.path import join

import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp
from tqdm import tqdm

from gsmmutils import MyModel, DATA_PATH
from gsmmutils.dynamic.cobrapy.cobrapy_dfba import CobrapyDFBA


def main():
    ts = np.linspace(0, 15, 10)  # Desired integration resolution and interval
    y0 = [0.13, 0.15]
    model = MyModel(join(DATA_PATH, 'models/model_dfba.xml'), "e_ActiveBiomass__cytop")
    cobrapy_dfba = CobrapyDFBA(model)
    with tqdm() as pbar:
        cobrapy_dfba.pbar = pbar
        sol = solve_ivp(
            fun=cobrapy_dfba.dynamic_system,
            # events=[cobrapy_dfba.infeasible_event],
            t_span=(ts.min(), ts.max()),
            y0=y0,
            t_eval=ts,
            rtol=1e-4,
            atol=1e-8,
            method='BDF'
        )
    print(sol)
    ax = plt.subplot(111)
    ax.plot(sol.t, sol.y.T[:, 0])
    ax2 = plt.twinx(ax)
    ax2.plot(sol.t, sol.y.T[:, 1], color='r')

    ax.set_ylabel('Biomass', color='b')
    ax2.set_ylabel('Phosphate', color='r')

    plt.show()


if "__main__" == __name__:
    main()