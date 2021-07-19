""" This script can train different machine learning models performs grid search and evaluation of models

"""


# Bash wrapper script to run machine_learning_tool.py in ml_docker_container
#
#     declare -a classifiers=("LinearSVC" "SVC" "RandomForest" "MultiLayerPerceptron" )
#
#     for val in ${classifiers[@]}; do \
#         echo Starting to train $val classifier \
#         python app/machine_learning_tool.py \
#         -i  <INPATH.TXT> \
#         -m <METADATA.CSV> \
#         -a $val \
#         -o <DEST-DIR>  \
#         -t <TITLE> \
#     done

# imports logging
import click
import logging
from time import time
import os

# imports data handling
import pandas as pd
import numpy as np
import fileUtils.metadata_handling as mh

# imports scikit-learn
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV
from sklearn.model_selection import cross_val_score, cross_validate
from sklearn.metrics import plot_confusion_matrix
from sklearn.metrics import plot_roc_curve
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit, StratifiedKFold

# scores
from sklearn.metrics import accuracy_score, precision_score, recall_score, balanced_accuracy_score
from sklearn.metrics import classification_report, make_scorer, matthews_corrcoef

# classifiers
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC, LinearSVC
from sklearn.ensemble import RandomForestClassifier

# feature selection
from sklearn.feature_selection import VarianceThreshold

# plotting
import umap
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

# Create logger
logger = logging.getLogger('Machine Learning Tool')
# Create console handler
ch = logging.StreamHandler()
# Create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
# add ch to logger
logger.addHandler(ch)
logger.setLevel(logging.INFO)


@click.command()
@click.option('-i', '--inpath', prompt='path to gene_count table',
              help='path to gene_count table', required=True)
@click.option('-m', '--metadata', prompt='path to metadata table',
              help='Table with metadata in csv format', required=True)
@click.option('-a', '--algorithm', prompt='please enter a machine learning algorithm',
              help='please enter a machine learning algorithm', required=True,
              type=click.Choice(['LinearSVC', 'SVC', 'RandomForest', 'MultiLayerPerceptron'], case_sensitive=True))
@click.option('-k', '--kegg', help='provide path to file with KEGG genes of interest',
              required=False, default=None)
@click.option('-o', '--outpath', prompt='path for outputs',
              help='Path for output tables and plots. Extension is filled automatically.')
@click.option('-c', '--cores', help='please enter number of cpu cores to be used for search', default=2)
@click.option('-t', '--title', help='please enter a dataset title/tissue title')
def main(inpath, metadata, algorithm, kegg, outpath, title, cores):
    start_time = time()

    # Make sure output directory exists
    new_outpath = outpath + 'ml_results/'
    if not os.path.exists(new_outpath):
        os.mkdir(new_outpath)

    logger.info('Parse gene_expression table...')
    if kegg:
        logger.info('Applying KEGG filter...')
        kegg_genes_pathways = pd.read_csv(kegg, sep="\t")
        kegg_ensembl = kegg_genes_pathways.iloc[:, 0].tolist()
        kegg_name = kegg_genes_pathways.iloc[:, 1].tolist()
        id2name_dict = dict(zip(kegg_ensembl, kegg_name))
        expr_df, sample_names, gene_ids, feature_names = read_dataset(inpath, ensembl_filter=kegg_ensembl,
                                                                      id2name=id2name_dict)
    # No filtering for KEGG genes
    else:
        expr_df, sample_names, gene_ids, feature_names = read_dataset(inpath, ensembl_filter=None,
                                                                      id2name=None)

    logger.info('Parse metadata...')
    meta_dict, target_names, target, annotations = mh.read_metadata(metadata, sample_names)

    logger.info('Split training and testing data...')
    X_train, X_test, y_train, y_test = train_test_split(expr_df, target, stratify=target, test_size=0.2,
                                                        random_state=99)
    logger.info('Feature scaling and filtering out 0-variance genes...')
    X_train_f, X_test_f, f_names_red, gene_ids_red = feature_reduction(X_train, X_test, feature_names)

    # Define models
    # pass randomState Instance to RandomForest upon instantiation (following best practices sklearn)
    models = {'LinearSVC': LinearSVC(random_state=37), 'SVC': SVC(random_state=42),
              'RandomForest': RandomForestClassifier(random_state=np.random.RandomState(0)),
              'MultiLayerPerceptron': MLPClassifier(random_state=87)}

    logger.info('Performing grid search...')
    grid_search(model=algorithm, model_dict=models, X_train=X_train_f, y_train=y_train, X_test=X_test_f, y_test=y_test,
                gene_ids=gene_ids_red, feature_names=f_names_red, splits=5, n_jobs=cores, refitting=True,
                outpath=new_outpath, dataset_title=title)

    end_time = time()
    logger.info('Process finished in ' + str(round(end_time - start_time, 2)) + "sec")


def read_dataset(path, ensembl_filter, id2name):
    """ Function to parse a dataset

    :param path: to input gene-expression table (GeneID, GeneName, Sample1, Sample2, ...)
    :param ensembl_filter: input list of ensembl_gene_ids
    :param names: input list of gene_names
    :return: pandas.Dataframe
    :return: sample_names, list
    :return: gene_ids = features, list
    :return: gene_names = feature_names, list
    """

    data = pd.read_csv(path, sep="\t")
    df = data.T
    feature_names = df.iloc[1].tolist()
    df = df.drop(["GeneName"], axis=0)
    df.columns = df.iloc[0]
    df = df.drop(["GeneID"], axis=0)
    df.columns.name = None
    sample_names = df.index.tolist()
    gene_ids = df.columns.tolist()
    df_final = df.fillna(0)

    if ensembl_filter is not None:
        to_select = df.columns[df.columns.isin(ensembl_filter)]
        filtered_df = df[to_select].copy()
        copy = filtered_df.copy()
        df_final = copy.fillna(0)
        gene_ids = df_final.columns.tolist()
        feature_names = [v for e in gene_ids for k, v in id2name.items() if e == k]
        # print(len(gene_ids), len(feature_names), len(df_final.columns))

    return df_final, sample_names, gene_ids, feature_names


def feature_reduction(X_train, X_test, featureNames):
    """ Function to filter out 0-variance genes

    :param X_train: numpy.array of training data
    :param X_test: numpy.array of testing data
    :param featureNames: array with gene feature names
    :return: filtered and scaled X_train, numpy array
    :return: filtered and scaled X_test, numpy array
    :return: filtered feature_names, list
    :return: filtered gene_ids, list
    """

    scaler = MinMaxScaler()
    var_filter = VarianceThreshold(threshold=0)

    # get gene_ids for reducing after variance threshold
    gene_ids = np.array(X_train.columns.tolist())

    X_train_sc = scaler.fit_transform(X_train)
    X_test_sc = scaler.transform(X_test)

    X_train_filtered = var_filter.fit_transform(X_train_sc)
    X_test_filtered = var_filter.transform(X_test_sc)

    # get names of reduced feature_names
    mask = var_filter.get_support()
    feature_arr = np.array(featureNames)
    feature_names_reduced = feature_arr[mask].tolist()
    gene_ids_red = gene_ids[mask].tolist()

    # print INFO
    print("{}: {}".format("Number of training instances ", X_train_filtered.shape[0]))
    print("{}: {}".format("Number of test instances ", X_test_filtered.shape[0]))
    print("{}: {}".format("Number of features after VarianceFilter ", len(feature_names_reduced)))
    print("{}: {}".format(len(featureNames) - len(feature_names_reduced), "features were removed "))

    return X_train_filtered, X_test_filtered, feature_names_reduced, gene_ids_red


def max_features_arr(feature_names, max_f_arr):
    """ Function for obtaining max_features as float array for RandomForest

    :param feature_names: array of gene feature names
    :param max_f_arr: array of max_features to be converted into float array for grid search
    :return: float array of max_features
    """

    n_features_ = len(feature_names)
    float_arr = [i / n_features_ for i in max_f_arr]
    return float_arr


def search_space(model):
    """Function to define the search space for GridSearch or RandomGridSearch

    :param model: one of LinearSVC, SVC, RandomForest, MultiLayerPerceptron
    :return: search space
    """

    space = {}

    if model == 'LinearSVC':
        # combinations = 26 * 7 = 182
        space = {'C': [0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 1, 2, 4, 6, 8, 10, 12, 16, 20, 24, 28, 32, 36, 40, 44,
                       48, 52, 56, 60, 64],
                 'max_iter': [1000, 2000, 5000, 10000, 20000, 50000, 100000]}

    if model == 'SVC':
        # combinations = 26 * 6 * 3 * 7 = 3276
        space = {'C': [0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 1, 2, 4, 6, 8, 10, 12, 16, 20, 24, 28, 32, 36, 40, 44,
                       48, 52, 56, 60, 64],
                 'gamma': [0.0001, 0.001, 0.01, 0.1, 1.0, 10.0],
                 'kernel': ['sigmoid', 'poly', 'rbf'],
                 'max_iter': [1000, 2000, 5000, 10000, 20000, 50000, 100000]}
    if model == 'RandomForest':
        # combinations = 9 * 2 * 32 * 3 * 3 * 2 = 10368
        space = {'n_estimators': [5, 10, 20, 50, 100, 200, 400, 600, 1000],
                 'max_features': ['auto', 'sqrt'],
                 'max_depth': np.linspace(1, 32, 32, endpoint=True),
                 'min_samples_leaf': [3, 4, 5],
                 'min_samples_split': [8, 10, 12],
                 'bootstrap': [True, False]
                 }
    if model == 'MultiLayerPerceptron':
        # combinations = 4 * 2 * 2 * 4 * 5 * 3 = 960
        space = {'hidden_layer_sizes': [(500, 100), (100, 50), (50, 20), (30, 10)],
                 'activation': ['tanh', 'relu'],
                 'solver': ['sgd', 'adam'],
                 'alpha': [0.0001, 0.001, 0.01, 0.1],
                 'max_iter': [1000, 5000, 10000, 20000, 50000, 100000],
                 'learning_rate': ['constant', 'adaptive']}

    return space


def grid_search(model, model_dict, X_train, y_train, X_test, y_test, gene_ids, feature_names, splits, n_jobs, refitting,
                outpath, dataset_title):
    """Function to perform grid search and print evaluation performance

    :param model: LinearSVC, SVC, RandomForest, MultiLayerPerceptron
    :param model_dict: dictionary which instantiates models
    :param X_train: scaled training data, numpy array
    :param y_train: list of target values
    :param X_test:  scaled test data, numpy array
    :param y_test:  list of target test values
    :param gene_ids: list of gene ids
    :param feature_names: list of features
    :param splits: number of splits for cross-validation
    :param n_jobs: number of cpu-cores to be used
    :param refitting: default is True
    :param outpath: for plots and tables
    :param dataset_title: i.e. tissue_dataset
    :return: None
    """

    # set display options
    pd.set_option("display.max_rows", 20, "display.max_columns", 15)

    # construct output path
    new_outpath = outpath + model + '_' + dataset_title

    # Select classifier from dict
    clf = model_dict[model]

    # Select grid
    grid = search_space(model)

    print("Performing random search for: ", clf)
    t0 = time()

    # Define scoring metrics
    scoring = {'mcc': make_scorer(matthews_corrcoef), 'bal_acc': 'balanced_accuracy', 'rec': 'recall',
               'prec': 'precision'}
    Matthews = make_scorer(matthews_corrcoef, greater_is_better=True)

    # Instantiate stratified-shuffled 5x2 cross validation: shuffling important when class labels contiguously in data
    cv = StratifiedKFold(n_splits=2, shuffle=True, random_state=65)
    # Define GridSearch
    search = RandomizedSearchCV(clf, grid, cv=cv, return_train_score=True, scoring=Matthews, n_jobs=n_jobs,
                                refit=refitting)  # when refit = True, fitting once to get  best_model is sufficient

    # Performing non-nested grid search and hyper-parameter tuning
    search.fit(X_train, y_train)
    print("done in %0.3fs" % (time() - t0))
    print("Best params: ", search.best_params_)
    print("Best estimator: ", search.best_estimator_)

    best_model = search.best_estimator_
    # Cross_val_scores of trees
    cv_results = search.cv_results_
    cv_df = pd.DataFrame.from_dict(cv_results)
    print('Cross-Validation GridSearch ' + model + '\n', cv_df)
    # print(cv_results)
    print_mean_std_scores(cv_results)

    # Outer CV loop (5x2 cross validation) with different scores
    nested_score = cross_validate(best_model, X_train, y_train, cv=splits, scoring=scoring)
    print("Outer cross validation scores on train data: ")
    print_mean_std_scores(nested_score)

    if model == 'RandomForest':
        # Extract feature importance of forest
        forest_feature_importance = best_model.feature_importances_
        df_forest_feature_importance = pd.DataFrame(sorted(zip(forest_feature_importance, gene_ids, feature_names),
                                                           reverse=True),
                                                    columns=['Feature_Importance', 'GeneID', 'GeneName'])
        plot_bar(df_forest_feature_importance, new_outpath)

        print("20 most important features RandomForest \n", df_forest_feature_importance.head(n=20))
        # Save feature_importance to file
        df_forest_feature_importance.to_csv(new_outpath + '_feature_importance.csv', index=False)

    # Predict on test data
    y_pred = best_model.predict(X_test)

    # Evaluate predictions
    print(model, ' : ')
    print(classification_report(y_test, y_pred, target_names=['normal', 'tumor']))
    print(model + " MCC on test data %.3f" % (matthews_corrcoef(y_test, y_pred)))
    print("Balanced Accuracy on test data %.3f" % (balanced_accuracy_score(y_test, y_pred)))

    # Evaluate model
    fig = plot_confusion_matrix(best_model, X_test, y_test, display_labels=('control', 'tumor'), cmap=plt.cm.Blues)
    fig.ax_.set_title('Confusion Matrix - ' + model)
    plt.savefig(new_outpath + '_confusion_matrix.png')

    roc_fig = plot_roc_curve(best_model, X_test, y_test, name=model)
    plt.savefig(new_outpath + '_roc_auc.png')


def print_mean_std_scores(scoring_dict):
    """Print cross_validation_results

    :param scoring_dict: from grid_search, cross_val_score, cross_validate
    :return: None
    """

    for k, v in scoring_dict.items():
        # parameters cannot be used for computing mean and std
        if k.startswith('param'):
            continue
        else:
            print("%s: mean=%.3f, (+/- std=%.3f)" % (k, v.mean(), v.std()))


def plot_bar(df, outpath):
    """Plot 10 most important features

    :param df: pandas Dataframe of feature importances
    :return: None
    """

    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

    # set the style of the axes and the text color
    plt.rcParams['axes.edgecolor'] = '#333F4B'
    plt.rcParams['axes.linewidth'] = 0.8
    plt.rcParams['xtick.color'] = '#333F4B'
    plt.rcParams['ytick.color'] = '#333F4B'
    plt.rcParams['text.color'] = '#333F4B'

    names = df.iloc[:, 2].tolist()[:10]
    percentages = df.iloc[:, 0].tolist()[:10]

    fig, ax = plt.subplots(figsize=(7, 5))
    plt.bar(names, percentages, width=0.5, color=colors['gray'], edgecolor=colors['black'])
    plt.xticks(rotation=45, ha='right')
    ax.set_ylabel('Feature Importance', fontsize=12)
    ax.set_xlabel('Genes', fontsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outpath + '_feature_importance_plot.png')


if __name__ == "__main__":
    main()
