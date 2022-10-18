
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def barplot(df, to_show = True, path = None):
    sns.set(rc = {'figure.figsize':(20,10)})
    df.plot(kind="bar", fontsize=12)
    plt.tight_layout(pad=2.5)
    plt.xlabel("condition", fontsize=14)
    plt.ylabel("growth rate $(d^{-1})$", fontsize=14)
    if to_show:
        plt.show()
    else:
        plt.savefig(path)