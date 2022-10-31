
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from scipy import stats
sns.set(rc={'figure.figsize': (20, 10)})


def barplot(dataframe, to_show = True, path = None):
    dataframe.plot(kind="bar", fontsize=12)
    plt.tight_layout(pad=2.5)
    plt.xlabel("condition", fontsize=14)
    plt.ylabel("growth rate $(d^{-1})$", fontsize=14)
    if to_show:
        plt.show()
    else:
        plt.savefig(path)


def boxplot(dataframe, x_cols=None, y_cols = None, to_show = True, path = None, x_labels = None, y_labels = None):
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


def surface(dataframe, x, y, z):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_trisurf(dataframe[y], dataframe[x], dataframe[z], cmap=plt.cm.jet, linewidth=0.2)
    plt.show()


def hist(dataset, columns, to_show = True, path = None, title = None, xlabel = None, ylabel = None):
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    sns.distplot(dataset[columns])
    if to_show:
        plt.show()
    else:
        plt.savefig(path)


def qqplot(model, to_show = True, path = None):
    res = model.resid
    fig = sm.qqplot(res, stats.t, fit=True, line="45")
    if to_show:
        plt.show()
    else:
        plt.savefig(path)


def heatmap(dataframe, to_show = True, path = None):
    sns.heatmap(dataframe, cmap="YlGnBu", z_score=0)
    if to_show:
        plt.show()
    else:
        plt.savefig(path)

def clustermap(dataframe, to_show = True, path = None):
    g = sns.clustermap(dataframe, cmap="Blues", z_score=0, center=0)
    if to_show:
        plt.show()
    else:
        plt.savefig(path)
    return g