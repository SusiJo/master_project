#!/usr/bin/python3

# The purpose of this script is to
# 1) Transform the merged TPM data such that is can by used by ML: row = sample, col = gene (feature)
# 2) Scale the data with MinMaxScaler from scikit-learn on genes = columns
# 3) Create a dictionary like structure for storage target_names = label [tumor, control]
# {data: array[], target: array[], target_names: array[], feature_names: array[]}


# imports
import os
import sys
import click
import numpy as np
import pandas as pd
from fileUtils.file_handling import read_tpm
import plotly.express as px
from time import time
from sklearn import preprocessing
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap # umap-learn
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from collections import OrderedDict


@click.command()
@click.option('-i', '--inpath', prompt='merged TPM table .txt',
              help='Path to file with TPM values',
              required=True)
@click.option('-m', '--metadata', prompt='metadata to TPM file',
              help='Path to file with metadata',
              required=True)
@click.option('-o', '--outpath', prompt='path to image folder', required=True)
@click.option('--pca/--no-pca', help='Perform PCA analysis', default=False)
@click.option('--tsne/--no-tsne', help='Perform TSNE analysis', default=False)
@click.option('--umap/--no-umap', help='Perform UMAP analysis', default=False)
@click.option('--comparison/--no-comparison', help='Perform Comparison analysis', default=False)
@click.option('--silhouette/--no-silhouette', help='Perform Silhouette analysis', default=False)
def main(inpath, outpath, metadata, pca, tsne, umap, comparison, silhouette):
    print("Reading data...")
    sample_ids, gene_ids, feature_names, data = read_tpm(inpath)
    # print("sample ids", sample_ids)
    # print("data", data)

    # convert data to numpy array for scaling and transpose to scale on genes
    # shape = (features, samples)
    new_data = np.array(data, dtype=float)

    # Transpose data for Dimensionality reduction
    # shape = (samples, features)
    print("Transposing data...")
    dataT = np.transpose(new_data)

    # scaling on the "columns" = per gene feature such that gene expr. values between (-1,1)
    print("Scaling data...")
    scaled_data = scaling(dataT)

    print("Reading metadata...")
    meta, target_names, target, annotations, project, color_list = read_metadata(metadata, sample_ids)

    plotly_comparison(scaled_data, project, annotations, target, outpath)

    # Plotting
    if pca:
        visualization_plots(scaled_data, annotations, outpath, 'pca', 'Pancreas - TPM from featureCounts ComBatSeq  outlier removed')
    if tsne:
        visualization_plots(scaled_data, annotations, outpath, 'tsne', 'Pancreas - TPM from featureCounts ComBatSeq  outlier removed')
    if umap:
        visualization_plots(scaled_data, annotations, outpath, 'umap', 'Pancreas - TPM from featureCounts ComBatSeq outlier removed')

    if comparison:
        comparison_dim_reduction(scaled_data, color_list, outpath)

    if silhouette:
        silhouette_plot(scaled_data, target)
        print("Silhouette Coefficient: ", silhouette_score(scaled_data, target))
        print(silhouette_samples(scaled_data, target))


def read_metadata(metadata, sample_ids):
    """Read in metadata file csv format
    FileID, SampleType, CaseID, Project

    :return: meta_dict
    :return: target_names: dict
    :return: target: array
    :return: annotation: pd.DataFrame
    :return: project_arr
    :return: color_list
    """
    meta_dict = {}
    with open(metadata, 'r') as file:
        for line in file.readlines()[1:]:
            fields = line.strip().split(",")
            file_id = fields[0].strip()
            case = fields[1].strip()
            sample_type = fields[2].strip()
            project = fields[3].strip()
            if file_id not in meta_dict:
                meta_dict[file_id] = [sample_type, case, project]

    target_names = {}
    for i, e in enumerate(sample_ids):
        if e not in target_names.keys():
            target_names[e] = meta_dict[e][0]

    target = [0 if value == 'normal' else 1 for key, value in target_names.items()]

    # Create Dataframe object for annotations for plots
    # sort the metadata
    sorted_meta = OrderedDict([(el, meta_dict[el]) for el in target_names])
    # print("sorted meta ", sorted_meta)
    tmp = {0 + i: [k, sorted_meta.get(k)] for i, k in enumerate(sorted_meta)}

    annos = pd.DataFrame.from_dict(tmp, orient='index')
    # expand list with metadata into own series
    tags = annos.iloc[:, 1].apply(pd.Series)
    tags.rename({0: 'Condition', 1: 'CaseID', 2: 'Project'}, axis=1, inplace=True)
    # print("tags ", tags)
    # concat dataframes
    annotations = pd.concat([annos[:], tags[:]], axis=1)
    annotations.drop([1], axis=1, inplace=True)
    annotations.rename({0: 'ID'}, axis=1, inplace=True)

    project_list = annotations.iloc[:, 3]
    # print(project_list)
    project_set = set(project_list)
    # print(project_set)
    project_dict = {v: 0 + k for k, v in enumerate(project_set)}
    # print(project_dict)
    project_arr = [project_dict[k] for k in project_list]
    # print(project_arr)
    c = ['b', 'g', 'r', 'p']
    color_tups = list(zip(project_set, c))
    color_dict = {e[0]: e[1] for e in color_tups}
    print(color_dict)
    color_list = [color_dict[k] for k in project_list]

    # print(annotations)
    # print(color_list)

    return meta_dict, target_names, target, annotations, project_arr, color_list


def scaling(data):
    """Scaling with MinMax to range (-1,1)

    :param data: numpy array
    :return: scaled_data
    """
    scaler = preprocessing.MinMaxScaler()
    # scales in the range of -1 to 1
    scaled_data = scaler.fit_transform(data, (-1, 1))
    return scaled_data


def visualization_plots(data, target, outpath, method=str, title_dataset=str):
    """Creates interactive plotly graphs with reduced dimensions with tooltip displaying metadata

    :param data: numpy array
    :param target: list
    :param outpath: to folder for images .png, .html
    :param method: one of [pca, tsne, umap]
    :param title_dataset: i.e. Pancreas ComBat corrected
    """

    if method == 'pca':
        print("Performing pca for plotting...")
        pca = PCA(n_components=2)
        X_pca = pca.fit_transform(data)
        d1 = X_pca[:, 0]
        d2 = X_pca[:, 1]
        evr1 = pca.explained_variance_ratio_[0]
        evr2 = pca.explained_variance_ratio_[1]

    if method == 'tsne':
        print("Performing tsne for plotting...")
        # Dimension reduction
        # consider setting init='pca'
        tsne = TSNE(n_components=2, random_state=0)
        X_embedded = tsne.fit_transform(data)
        d1 = X_embedded[:, 0]
        d2 = X_embedded[:, 1]

    if method == 'umap':
        print("Performing umap for plotting...")
        # Dimension reduction
        # ensure reproducibility of embedding by setting random state
        reducer = umap.UMAP(random_state=42)
        embedding = reducer.fit_transform(data)
        d1 = embedding[:, 0]
        d2 = embedding[:, 1]

    # Plotting
    fig = px.scatter(x=d1, y=d2, color=target.Project, hover_name=target.ID, hover_data={'condition': target.Condition, 'case_id': target.CaseID})    # color project
    # fig = px.scatter(x=d1, y=d2, color=target.Condition, hover_name=target.ID, hover_data={'project': target.Project, 'case_id': target.CaseID})        # color condition
    if method == 'pca':
        fig.update_layout(title=str(method).upper() + ' - ' + title_dataset + ' dataset',
                          xaxis_title="PC1 - explained variance: " + str(round(evr1*100, 2)),
                          yaxis_title="PC2 - explained variance: " + str(round(evr2*100, 2)),
                          legend_title="Project")   # Condition
    else:
        fig.update_layout(title=str(method).upper() + ' - ' + title_dataset + ' dataset',
                          xaxis_title="Dim1", yaxis_title="Dim2", legend_title="Project")     # Project

    fig.update_traces(marker=go.scatter.Marker(size=20))
    fig.show()
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    # fig.write_image(outpath + str(method) + "_tpm_from_featureCounts_condition_plot.png")
    fig.write_html(outpath + str(method) + "_tpm_combatseq_no_outlier_plot.html")


def comparison_dim_reduction(X, target, outpath):
    """Main code from scikit-learn to create comparative plots with different reduction methods

    :param X: data
    :param target: list
    :param outpath: to store image .png
    """
    # Create figure
    fig = plt.figure(figsize=(15, 5))
    fig.suptitle("Dimensionality reduction with ", fontsize=14)
    n_components = 2

    # Set-up manifold methods
    methods = OrderedDict()
    methods['PCA'] = PCA(n_components=n_components, random_state=0)
    # init='pca'
    methods['t-SNE'] = TSNE(n_components=n_components, random_state=0)
    methods['UMAP'] = umap.UMAP(random_state=42)

    # Plot results
    for i, (label, method) in enumerate(methods.items()):
        t0 = time()
        Y = method.fit_transform(X)
        t1 = time()
        print("%s: %.2g sec" % (label, t1 - t0))
        ax = fig.add_subplot(1, 3, 1 + i)
        ax.scatter(Y[:, 0], Y[:, 1], c=target) # , cmap=plt.cm.Spectral), cmap=plt.cm.Set2
        ax.set_title("%s (%.2g sec)" % (label, t1 - t0))
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.yaxis.set_major_formatter(NullFormatter())
        ax.axis('tight')

    # plt.legend()
    # plt.show()
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    plt.savefig(outpath + "comparison_plot.png", format='png')


def plotly_comparison(data, target, annos, conditions, outpath):
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(data)
    d1 = X_pca[:, 0]
    d2 = X_pca[:, 1]
    evr1 = pca.explained_variance_ratio_[0]
    evr2 = pca.explained_variance_ratio_[1]

    tsne = TSNE(n_components=2, random_state=0)
    X_embedded = tsne.fit_transform(data)
    t1 = X_embedded[:, 0]
    t2 = X_embedded[:, 1]

    reducer = umap.UMAP(random_state=42)
    embedding = reducer.fit_transform(data)
    u1 = embedding[:, 0]
    u2 = embedding[:, 1]

    fig = make_subplots(rows=2, cols=3, subplot_titles=("PCA", "t-SNE", "UMAP"))
    fig.add_trace(go.Scatter(x=d1, y=d2, mode='markers', marker=dict(color=target)), row=1, col=1)
    fig.add_trace(go.Scatter(x=t1, y=t2, mode='markers', marker=dict(color=target)), row=1, col=2)
    fig.add_trace(go.Scatter(x=u1, y=u2, mode='markers', marker=dict(color=target)), row=1, col=3)
    fig.update_traces(hovertext='ID: ' + annos.ID + '<br>' + 'Project: ' + annos.Project + '<br>' + 'Condition: ' + annos.Condition + '<br>' + 'CaseID: ' + annos.CaseID)
    # fig.update_traces(marker_color=target)
    fig.add_trace(go.Scatter(x=d1, y=d2, mode='markers', marker=dict(color=conditions)), row=2, col=1)
    fig.add_trace(go.Scatter(x=t1, y=t2, mode='markers', marker=dict(color=conditions)), row=2, col=2)
    fig.add_trace(go.Scatter(x=u1, y=u2, mode='markers', marker=dict(color=conditions)), row=2, col=3)
    fig.update_traces(hovertext='ID: ' + annos.ID + '<br>' + 'Project: ' + annos.Project + '<br>' + 'Condition: ' + annos.Condition + '<br>' + 'CaseID: ' + annos.CaseID)
    # fig.update_traces(marker_color=conditions)
    fig.update_layout(showlegend=False, width=2000, height=1200, title_text='Dimensionality Reduction Comparison ComBat corrected')
    fig.show()

    if not os.path.exists(outpath):
        os.mkdir(outpath)
    fig.write_html(outpath + "comparison_combat.html")


def silhouette_plot(data, target):
    """Determine good number of clusters

    :param data: numpy array
    :param target: list
    """
    fig, ax = plt.subplots()
    ax.set_title("The silhouette plot")
    ax.set_xlabel("The silhouette coefficient values")
    ax.set_ylabel("Cluster label")

    # The subplot is the silhouette plot
    ax.set_xlim([-0.1, 1])

    lower = 10
    clusters = ['normal', 'tumor']
    target = np.array(target)

    # UMAP
    reducer = umap.UMAP(random_state=10)
    embedding = reducer.fit_transform(data)

    # tSNE
    # tsne = TSNE(n_components=2, random_state=0)
    # embedding = tsne.fit_transform(data)

    for i, condition in enumerate(clusters):

        # Compute the silhouette score average
        silhouette_avg = silhouette_score(embedding, target)
        # Compute the silhouette scores for each sample
        sample_silhouette_values = silhouette_samples(embedding, target)

        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = sample_silhouette_values[target == i]
        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]

        upper = lower + size_cluster_i

        color = cm.nipy_spectral(float(i) / 2)
        # plotting the silhouette clusters
        ax.fill_betweenx(np.arange(lower, upper),
                         0, ith_cluster_silhouette_values,
                         facecolor=color, edgecolor=color, alpha=0.7)

        # Label the silhouette plots with their cluster numbers at the middle
        ax.text(-0.05, lower + 0.5 * size_cluster_i, condition)

        # shift by 10 to plot space between silhouettes
        lower = upper + 10  # 10 for the 0 samples

    # The vertical line for average silhouette score of all the values
    ax.axvline(x=silhouette_avg, color="red", linestyle="--")
    ax.set_yticks([])  # Clear the yaxis labels / ticks
    ax.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

    fig.show()


if __name__ == "__main__":
    sys.exit(main())
