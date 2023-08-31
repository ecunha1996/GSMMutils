from typing import Union, List, Any, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats

sns.set(rc={'figure.figsize': (15, 6)})


def barplot(dataframe, to_show=True, path=None):
    dataframe.plot(kind="bar", fontsize=12)
    plt.tight_layout(pad=2.5)
    plt.xlabel("condition", fontsize=14)
    plt.ylabel("growth rate $(d^{-1})$", fontsize=14)
    if to_show:
        plt.show()
    else:
        plt.savefig(path)


def lineplot(x, y, to_show=True, **kwargs):
    graph = (sns.lineplot(x=x, y=y))
    graph.set(**kwargs)
    if to_show:
        plt.show()
    return graph

def boxplot(dataframe, x_cols=None, y_cols=None, to_show=True, path=None, x_labels=None, y_labels=None):
    for y in y_cols:
        fig, axs = plt.subplots(ncols=4)
        fig.tight_layout(pad=3.0)
        index = 0
        for x in x_cols:
            sns.boxplot(data=dataframe, x=x, y=y, ax=axs[index]).set(xlabel=x_labels[x], ylabel=y_labels[y])
            index += 1
    if to_show:
        plt.show()
    else:
        plt.savefig(path)


def surface(dataframe, x, y, z, **kwargs):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_trisurf(dataframe[y], dataframe[x], dataframe[z], linewidth=0.2, **kwargs)
    plt.show()


def hist(dataset, columns, to_show=True, path=None, title=None, xlabel=None, ylabel=None):
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    print(dataset[columns])
    sns.distplot(dataset[columns])
    if to_show:
        plt.show()
    else:
        plt.savefig(path)


def qqplot(model, to_show=True, path=None):
    import statsmodels.api as sm
    res = model.resid
    sm.qqplot(res, stats.t, fit=True, line="45")
    if to_show:
        plt.show()
    else:
        plt.savefig(path)


def heatmap(dataframe, to_show=True, path=None):
    sns.heatmap(dataframe, cmap="YlGnBu")
    if to_show:
        plt.show()
    else:
        plt.savefig(path)


def clustermap(dataframe, to_show=True, path=None, title="Title", **kwargs):
    g = sns.clustermap(dataframe, cmap="Blues", z_score=0, center=0, **kwargs).fig.suptitle(title)
    if to_show:
        plt.show()
    else:
        plt.savefig(path)
    return g


def plot_two_axis(df, secondary=None):
    if type(df) == pd.Series:
        df = df.to_frame()
    if secondary is None:
        secondary = []
    _, ax = plt.subplots(figsize=(20, 10))
    colors = {'C00009__extr': 'red', 'C00244__extr': 'green', 'C02094__chlo': "black"}
    for column in df.columns:
        if column not in secondary:
            ax.plot(df.index, df[column], label=column, color=colors[column])
        else:
            ax2 = ax.twinx()
            ax2.plot(df.index, df[column], label=column)
    plt.show()


def plot_concentrations(data: pd.DataFrame(), y: Union[str, list] = None, to_show: bool = False,
                        experimental: Union[List[List], List[Tuple[List, Any]]] = None, filename=None,
                        x_label=None, y_label=None, title=None, secondary_axis: Union[str, list] = None,
                        secondary_y_label=None, experimental_label=None):
    try:
        if not x_label:
            x_label = 'time (d)'
        if not y_label:
            y_label = 'concentration (g/L)'
        if type(y) == str:
            y = [y]
        if type(secondary_axis) == str:
            secondary_axis = [secondary_axis]
        plt.style.context('Solarize_Light2')
        fig = plt.figure()
        # plot all concentrations from the dataframe
        ax = plt.subplot(111)
        for i in y:
            ax.plot(data['time'], data[i])
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(title)
        legend = y
        if experimental:
            for index, exp in enumerate(experimental):
                ax.scatter(exp[0], exp[1])
                if not experimental_label:
                    legend += ["Experimental"]
                else:
                    legend += [experimental_label[index]]
        if secondary_axis:
            ax2 = ax.twinx()
            for i in secondary_axis:
                ax2.plot(data['time'], data[i], color='r')
            ax2.set_ylabel(secondary_y_label)
            legend += secondary_axis
            fig.add_subplot(ax2)
        fig.add_subplot(ax)
        fig.legend(legend, loc='upper right', bbox_to_anchor=(1.0, 1.0), bbox_transform=ax.transAxes)
        if to_show:
            fig.show()
        else:
            fig.savefig(filename)
        return fig
    except Exception as e:
        print(e)


def plot_trajectories(data: pd.DataFrame(), y: str = None, to_show: bool = False):
    fig = plt.figure()

    ax = plt.subplot(111)
    ax.plot(data['time'], data[y])

    ax2 = plt.twinx(ax)
    ax2.plot(data['time'], data['Flux'], color='r')

    ax.set_ylabel(y, color='b')
    ax2.set_ylabel('EX_e_Biomass__dra', color='r')
    ax2.scatter([0, 2, 4, 6, 8, 10, 12], [0.126, 0.221, 0.371, 0.458, 0.532, 0.576, 0.581])

    fig.add_subplot(ax)
    if to_show:
        fig.show()
    else:
        fig.savefig('concentrations.png')
    return fig


def generate_plot_for_data(filename: str, experimental_data: dict, simulation_data: dict, sd_data: dict, y_label: str):
    """
    Generates a plot for a given data
    Parameters
    ----------
    filename: str
        Filename for the plot
    experimental_data: dict
        Experimental data
    simulation_data: dict
        Simulation data
    sd_data: dict
        Standard deviation data
    y_label: str
        Label for the y axis

    Returns
    -------

    """
    plt.clf()
    ax = plt.subplot(111)
    sns.barplot(x=list(experimental_data.keys()), y=list(experimental_data.values()), errorbar='sd')
    if sd_data:
        plt.errorbar(x=list(experimental_data.keys()), y=list(experimental_data.values()), yerr=list(sd_data.values()),
                     fmt='none', color='black', capsize=4)
    ax.scatter(x=list(simulation_data.keys()), y=list(simulation_data.values()), zorder=2)
    plt.ylabel(y_label)
    plt.xlabel(r"Trial")
    plt.savefig(filename)
