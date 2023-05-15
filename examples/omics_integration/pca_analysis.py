from typing import Tuple

import pandas as pd
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler

from matplotlib import colors as mcolors
import prince

from ExpGSMM import DATA_PATH
from Tissue_specific_Reconstruction_Pipeline.Pipeline.utils.config_variables import DATASET

M_COLORS = list(mcolors.TABLEAU_COLORS.keys())
COLORS = M_COLORS + ['#1f6357', '#017374', '#0cb577', '#ff0789', '#afa88b']


def get_explained_variance_idx(component):
    if '1' in component:
        return 0

    elif '2' in component:
        return 1

    elif '3' in component:
        return 2

    return 3


def scaling(dataframe: pd.DataFrame,
            categorical: list,
            standard: bool = True,
            variance: bool = True) -> pd.DataFrame:
    mask = dataframe.columns.isin(categorical)
    cols = dataframe.columns[~mask]
    x = dataframe.loc[:, cols]
    y = dataframe.loc[:, categorical]

    if variance:
        scalier = VarianceThreshold()
        scaled = scalier.fit(x)
        x = x.iloc[:, scaled.get_support(indices=True)]

    if standard:
        scalier = StandardScaler()
        scaled = scalier.fit_transform(x.T)

        x = pd.DataFrame(scaled, columns=x.index, index=x.columns)
        x = x.T

    return pd.concat([x, y], axis=1)


def pca_analysis(dataframe: pd.DataFrame,
                 categorical: list,
                 components: int = 2) -> Tuple[pd.DataFrame, PCA]:
    mask = dataframe.columns.isin(categorical)
    cols = dataframe.columns[~mask]
    x = dataframe.loc[:, cols]
    y = dataframe.loc[:, categorical]

    pca = PCA(n_components=components)

    pc = pca.fit_transform(x)

    columns = [f'PC {i + 1}' for i in range(components)]

    df = pd.DataFrame(data=pc, index=dataframe.index, columns=columns)

    df = pd.concat([df, y], axis=1)
    return df, pca


def run_mca(n_components: int, data: pd.DataFrame, categorical: list):
    mask = data.columns.isin(categorical)
    cols = data.columns[~mask]

    x = data.loc[:, cols]
    y = data.loc[:, categorical]

    x = x.astype(str)

    mca = prince.MCA(n_components=n_components)

    mc = mca.fit(x)
    df_mca = mc.transform(x)

    df_mca.columns = [f'PC {i + 1}' for i in range(n_components)]
    df = pd.concat([df_mca, y], axis=1)

    explained_inertia = mca.explained_inertia_
    return df, explained_inertia


def plot_pca(workdir: str,
             dataframe: pd.DataFrame,
             pca: PCA,
             c1: str,
             c2: str,
             title: str,
             factor: str):
    explained_variance_1 = get_explained_variance_idx(c1)
    explained_variance_2 = get_explained_variance_idx(c2)

    fig = plt.figure(figsize=(8, 8))

    ax = fig.add_subplot(1, 1, 1)
    ax.set_title(title, fontsize=20)

    x_label = f'{c1} ({round(pca.explained_variance_ratio_[explained_variance_1] * 100, 2)} %)'
    y_label = f'{c2} ({round(pca.explained_variance_ratio_[explained_variance_2] * 100, 2)} %)'
    ax.set_xlabel(x_label, fontsize=15)
    ax.set_ylabel(y_label, fontsize=15)

    labels = set(dataframe.loc[:, factor])

    for label, color in zip(labels, COLORS):
        mask = dataframe.loc[:, factor] == label

        pc1 = dataframe.loc[mask, c1]
        pc2 = dataframe.loc[mask, c2]
        ax.scatter(pc1,
                   pc2,
                   c=color,
                   s=60)

        organisms_id = dataframe.loc[mask, 'Tissue_id']

        for pc1_pt, pc2_pt, annotation in zip(pc1, pc2, organisms_id):
            ax.annotate(annotation, (pc1_pt + 1, pc2_pt - 1))

    legend = ax.legend(labels, loc=(1.04, 0))
    ax.grid()
    fig.savefig(fname=workdir, bbox_extra_artists=(legend,), bbox_inches='tight', dpi=300)


def plot_mca(data: pd.DataFrame, explained_variance: list, c1: str, c2: str, title: str, workdir: str, factor: str):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title(title, fontsize=20)
    x_label = f'{c1} ({round(explained_variance[0] * 100, 2)} %)'
    y_label = f'{c2} ({round(explained_variance[1] * 100, 2)} %)'
    ax.set_xlabel(x_label, fontsize=15)
    ax.set_ylabel(y_label, fontsize=15)

    labels = set(data.loc[:, factor])

    for label, color in zip(labels, COLORS):
        mask = data.loc[:, factor] == label

        pc1 = data.loc[mask, c1]
        pc2 = data.loc[mask, c2]
        ax.scatter(pc1,
                   pc2,
                   c=color,
                   s=60)
        #
        # organisms_id = data.loc[mask, 'Condition']
        #
        # for pc1_pt, pc2_pt, annotation in zip(pc1, pc2, organisms_id):
        #     ax.annotate(annotation, (pc1_pt + 1, pc2_pt - 1))

    legend = ax.legend(labels, loc=(1.04, 0))
    ax.grid()
    fig.savefig(fname=workdir, bbox_extra_artists=(legend,), bbox_inches='tight', dpi=300)


if __name__ == '__main__':
    # contents = ['Flux', 'Reaction']
    contents = ['Reaction']

    for content in contents:
        df = pd.read_csv(f'{DATA_PATH}/omics/{DATASET}/Dsalina_light_fastcore_Local2_2_4_4.csv', index_col=[0])
        df = df*1
        conditions, threshold = [], []
        for i, index in enumerate(df.index):
            conditions.append(index.split('_')[1])
            threshold.append('.'.join(index.split('_')[-2:]))

        df['Condition'] = conditions
        df['Threshold'] = threshold

        categorical = ['Condition', 'Threshold']

        df = scaling(dataframe=df, categorical=categorical)
        pca_df, pca = pca_analysis(dataframe=df, categorical=categorical, components=3)
        mca_df, explained_inertia = run_mca(n_components=3, data=df, categorical=categorical)

        factor_list = categorical

        title = ''
        if content == 'Flux':
            for factor in factor_list:
                plot_pca(workdir=f'{DATA_PATH}/omics/{DATASET}/Desai_{content}_PCA_2_{factor}.png',
                         dataframe=pca_df,
                         pca=pca,
                         c1='PC 1',
                         c2='PC 2',
                         title=title,
                         factor=factor)

                plot_pca(workdir=f'{DATA_PATH}/omics/{DATASET}/Desai_{content}_PCA_2_{factor}.png',
                         dataframe=pca_df,
                         pca=pca,
                         c1='PC 1',
                         c2='PC 3',
                         title=title,
                         factor=factor)

        elif content == 'Reaction':
            for factor in factor_list:
                plot_mca(data=mca_df,
                         explained_variance=explained_inertia,
                         c1='PC 1',
                         c2='PC 2',
                         title=title,
                         workdir=f'{DATA_PATH}/omics/{DATASET}/Desai_{content}_MCA_1_{factor}.png',
                         factor=factor)

                plot_mca(data=mca_df,
                         explained_variance=explained_inertia,
                         c1='PC 1',
                         c2='PC 3',
                         title=title,
                         workdir=f'{DATA_PATH}/omics/{DATASET}/Desai_{content}_MCA_2_{factor}.png',
                         factor=factor)
