from collections import OrderedDict
import pandas as pd


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
    # Read metadata into dictionary
    meta_dict = {}
    with open(metadata, 'r') as file:
        for line in file.readlines()[1:]:
            fields = line.strip().split(",")
            file_id = fields[0].strip()
            case = fields[1].strip()
            sample_type = fields[2].strip()
            project = fields[3].strip()
            if file_id not in meta_dict:
                meta_dict[file_id] = [case, sample_type, project]

    # Map sample_ids to conditions to create condition array (str)
    target_names = {}
    for i, e in enumerate(sample_ids):
        if e not in target_names.keys():
            target_names[e] = meta_dict[e][1]

    # Array of conditions encoded as numbers
    target = [0 if value == 'normal' else 1 for key, value in target_names.items()]

    # Sort the metadata
    sorted_meta = OrderedDict([(el, meta_dict[el]) for el in target_names])
    tmp = {0 + i: [k, sorted_meta.get(k)] for i, k in enumerate(sorted_meta)}

    # Create annotation dataframe from sorted dictionary (ID, CaseID, Condition, Project)
    annos = pd.DataFrame.from_dict(tmp, orient='index')
    tags = annos.iloc[:, 1].apply(pd.Series)
    tags.rename({0: 'CaseID', 1: 'Condition', 2: 'Project'}, axis=1, inplace=True)
    annotations = pd.concat([annos[:], tags[:]], axis=1)
    annotations.drop([1], axis=1, inplace=True)
    annotations.rename({0: 'ID'}, axis=1, inplace=True)

    return meta_dict, target_names, target, annotations
