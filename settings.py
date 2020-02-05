from config import PERMUTATION_REPETITIONS

# ---------------------------------- #
# CPM tests
# ---------------------------------- #
from cdg.changedetection import Student_t_test, ChangePointMethod, EDivisive_R, MMDCPM, EnergyCPM
MUCPM = 'mucpm'
ECPM = 'ecpm'
DCPM = 'dcpm'
KCPM = 'kcpm'
EDIV = 'ediv'
LIST_CPM = [MUCPM, ECPM, DCPM, KCPM, EDIV]
ECPM_copy = EnergyCPM(repetitions=PERMUTATION_REPETITIONS)
ECPM_copy.name = 'DCPM'
CPM_TEST_INSTANCES = {MUCPM: ChangePointMethod(local_test=Student_t_test()),
                    ECPM: EnergyCPM(repetitions=PERMUTATION_REPETITIONS),
                    DCPM: ECPM_copy,
                    KCPM: MMDCPM(repetitions=PERMUTATION_REPETITIONS),
                    EDIV: EDivisive_R(repetitions=PERMUTATION_REPETITIONS)}


# ---------------------------------- #
# Graph measures
# ---------------------------------- #
from cdg.graph import GraphEditDistanceNX, ShortestPathKernel, FrobeniusGraphDistance #, WeisfeilerLehmanKernel
frob = FrobeniusGraphDistance()
gednx_neuc = GraphEditDistanceNX(node_cost='euclidean')
# wlkernel = WeisfeilerLehmanKernel(verbose=True, normalize=True)
spkernel = ShortestPathKernel(verbose=True, normalize=True)


# ---------------------------------- #
# Embedding methods
# ---------------------------------- #
from cdg.embedding import DissimilarityRepresentation
NO_EMBEDDING_USE_DISTANCE = lambda emb_dim: 'NO_EMBEDDING_USE_DISTANCE'
NO_EMBEDDING_USE_KERNEL = lambda emb_dim: 'NO_EMBEDDING_USE_KERNEL'
EMBEDDINGS = {}
for c in LIST_CPM:
    # embeddings[t] = MultiDimensionalScaling
    EMBEDDINGS[c] = DissimilarityRepresentation
EMBEDDINGS[DCPM] = NO_EMBEDDING_USE_DISTANCE
EMBEDDINGS[KCPM] = NO_EMBEDDING_USE_KERNEL


# ---------------------------------- #
# Datasets and associated metrics
# ---------------------------------- #
# Dataset paths
letter_folder = "./datasets/iam/Letter/HIGH"
aids_folder = "./datasets/iam/AIDS/data"
mutag_folder = "./datasets/iam/Mutagenicity/data"
delaunay_folder = './datasets/delaunay'
delaunay_var_folder = './datasets/delaunay_var'
kaggle_dog_folder = "./datasets/kaggle-seizure/dog"
kaggle_human_folder = "./datasets/kaggle-seizure/hum"
# Classes
delaunay_class1_list = [-1, 8, 9, 10, 11, 12]
delaunay_ratio_class1_list = [8, 9, 10, 11, 12]
kaggle_dog_list = range(1, 5)
kaggle_human_list = range(1, 9)
# dynamic_sbm_setups = ['noch5050', 'glob', 'prop', 'i1090', 'm5050', 'f5050']
# sbm_iid_setups = ['gllowden', 'glmidden', 'glhighden']
# sbm_iid_setups += ['gllowden_poi', 'glmidden_poi', 'glhighden_poi']

# List of data sets settings
class DatasetSetting(object):
    def __init__(self, id, mea, gformat, folder):
        self.id = id
        self.mea = mea
        self.gformat = gformat
        self.folder = folder

dataset_settings = {}
dataset_settings['DELGED']    = DatasetSetting(id='Delaunay',   mea=gednx_neuc, gformat='nx',     folder=delaunay_folder)
dataset_settings['DELKSP']    = DatasetSetting(id='Delaunay',   mea=spkernel,   gformat='grakel', folder=delaunay_folder)
dataset_settings['DVARGED']   = DatasetSetting(id='DelVar',     mea=gednx_neuc, gformat='nx',     folder=delaunay_var_folder)
dataset_settings['DVARKSP']   = DatasetSetting(id='DelVar',     mea=spkernel,   gformat='grakel', folder=delaunay_var_folder)
# dataset_settings['DCSBMFRAG'] = DatasetSetting(id='DCSBM_frag', mea=frob,       gformat='cdg',    folder='./datasets/dcsbm_frag')


# ---------------------------------- #
# List of experiement settings
# ---------------------------------- #
experiments_with_multiple_settings = ['del', 'del_var', 'del_ratio', 'dog', 'human'] #todo this shouldn't be necessary

class ExperimentSetting(object):
    def __init__(self, id, cls, ratio, folder, name, dynamic_setup=None):
        self.id = id
        self.cls = cls
        self.ratio = ratio
        self.folder = folder
        self.name = name
        self.dynamic_setup = dynamic_setup

experiment_settings = [
    ExperimentSetting(id='debug',     cls=[0],  ratio=1.,  folder=delaunay_folder,  name="DelDebug"),
    ExperimentSetting(id='del',       cls=[0],  ratio=1.,  folder=delaunay_folder,  name="Delaunay"),
    ExperimentSetting(id='del_ratio', cls=[0],  ratio=.67, folder=delaunay_folder,  name="DelaunayRatio"),
    ExperimentSetting(id='del_var',   cls=[30], ratio=1.,  folder=delaunay_var_folder, name="DelaunayVar"),
    ExperimentSetting(id='del_multi1', cls=[0, 8, 12],        ratio=1., folder=delaunay_folder, name="DelaunayMulti"),
    ExperimentSetting(id='del_multi2', cls=[0, 8, 9, 10, 11], ratio=1., folder=delaunay_folder, name="DelaunayMulti"),
    ExperimentSetting(id='letter_multi1', cls=['A', 'E', 'H'],           ratio=1., folder=letter_folder, name="Letter"),
    ExperimentSetting(id='letter_multi2', cls=['A', 'E', 'F', 'H', 'I'], ratio=1., folder=letter_folder, name="Letter"),
    ExperimentSetting(id='aids', cls=['i', 'a'],           ratio=1., folder=aids_folder, name="AIDS"),
    ExperimentSetting(id='mutag', cls=['nonmutagen', 'mutagen'], ratio=1., folder=mutag_folder, name="Mutagenicity"),
    ExperimentSetting(id='del_no', cls=[0], ratio=1., folder=letter_folder, name="DelaunayNoChange"),
    ExperimentSetting(id='aids_no', cls=['i'], ratio=1., folder=aids_folder, name="AIDSNoChange"),
    ExperimentSetting(id='mutag_no', cls=['nonmutagen'], ratio=1., folder=mutag_folder, name="MutagenicityNoChange")
    ]
for i in kaggle_dog_list:
    experiment_settings.append(ExperimentSetting(id='dog{}'.format(i), cls=['preictal', 'ictal'], ratio=1., folder="{}{}".format(kaggle_dog_folder, i), name="Dog{}".format(i)))
for i in kaggle_human_list:
    experiment_settings.append(ExperimentSetting(id='hum{}'.format(i), cls=['preictal', 'ictal'], ratio=1., folder="{}{}".format(kaggle_human_folder, i), name="Human{}".format(i)))



sbm_nc_list = []
ct = 0
for e in ['glob', 'prop']:
    for n in [40, 100]:
        for a in [ 1. - .1 * i for i in range(10)][::-1]:
            ct += 1
            sbm_nc_list.append('sbm{}_{}'.format(e, ct))
            experiment_settings.append(ExperimentSetting(id=sbm_nc_list[-1], cls=[None], ratio=1., folder=None,
                                                         name="SBM_nc_{}a.{}%v.{}".format(e, round(100*a), n),
                                                         dynamic_setup={'alpha_update': a,
                                                                        'distance': frob.get_measure_fun(),
                                                                        'setting': e,
                                                                        'no_nodes': n,
                                                                        'seq_lengths': [200, 0]}))

sbm_driftplot_list = []
ct = 0
for e in ['glob']:
    for n in [40]:
        for a in [ 1. - .1 * i for i in range(10)][::-1]:
            ct += 1
            sbm_driftplot_list.append('sbm{}_{}'.format(e, ct))
            experiment_settings.append(ExperimentSetting(id=sbm_driftplot_list[-1], cls=[None], ratio=1., folder=None,
                                                         name="SBM_nc_{}a.{}%v.{}".format(e, round(100*a), n),
                                                         dynamic_setup={'alpha_update': a,
                                                                        'distance': frob.get_measure_fun(),
                                                                        'setting': e,
                                                                        'no_nodes': n,
                                                                        'seq_lengths': [200, 0]}))

sbm_prop_glob_list = []
for e in ['glob', 'prop']:
    for n in [40, 100]:
        for a in [ 0.5, 0.75, 1.]:
            ct += 1
            sbm_prop_glob_list.append('sbm{}_{}'.format(e, ct))
            experiment_settings.append(ExperimentSetting(id=sbm_prop_glob_list[-1], cls=[None], ratio=1., folder=None,
                                                         name="SBM_{}a.{}%v.{}".format(e, round(100*a), n),
                                                         dynamic_setup={'alpha_update': a,
                                                                        'distance': frob.get_measure_fun(),
                                                                        'setting': e,
                                                                        'no_nodes': n,
                                                                        'seq_lengths': [200, 50]}))


sbm_imf_list = []
for e in ['i1090', 'm5050', 'f5050']:
    for n in [40]:
        for a in [ 0.5, 0.75, 1.]:
            ct += 1
            sbm_imf_list.append('sbm{}_{}'.format(e, ct))
            experiment_settings.append(ExperimentSetting(id=sbm_imf_list[-1], cls=[None], ratio=1., folder=None,
                                                         name="SBM_{}a.{}%v.{}".format(e, round(100*a), n),
                                                         dynamic_setup={'alpha_update': a,
                                                                        'distance': frob.get_measure_fun(),
                                                                        'setting': e,
                                                                        'no_nodes': n,
                                                                        'seq_lengths': [200, 50]}))
