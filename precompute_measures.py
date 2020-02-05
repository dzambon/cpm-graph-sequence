import argparse
import sys
import numpy as np
from datetime import datetime

from cdg.utils import logger as cdg_logger
from cdg.graph import DataSet, isdataset, has_prec_measure, KERNEL, DISTANCE

from settings import dataset_settings
available_experiments = dataset_settings.keys()

# Define the argument parser
parser = argparse.ArgumentParser(description='')
parser.add_argument('-e', '--experiment', type=str,
                    required=True,
                    choices=available_experiments,
                    help='experiment name')

parser.add_argument('--force', action='store_true', default=False,
                    help='if dataset is already present, this will regenerate from scratch.')
parser.add_argument('-f', '--folder', type=str, default=None,
                    help='folder of the dataset')
parser.add_argument('-j', '--noJobs', type=int, default=1,
                    help='number of jobs of joblib')


def main(argv, logger=None):
    if logger is None:
        logger = cdg_logger
        logger.set_stdout_level(logger.DEBUG)
        logger.set_filelog_level(logger.DEBUG)
        
    args = parser.parse_args(argv)

    # command to parse
    launchingCommand = ''
    for a in sys.argv[1:]:
        launchingCommand += ' ' + a
    logger.info("command: " + launchingCommand)

    # check experiment availability
    if not args.experiment in available_experiments:
        raise ValueError('experiment on {} not available'.format(args.experiment))

    # load dataset
    if args.folder is not None:
        dataset_folder = args.folder
    else:
        dataset_folder = dataset_settings[args.experiment].folder
    assert isdataset(dataset_folder), 'no data set at {}.'.format(dataset_folder)

    has_prec_mea, mea_type = has_prec_measure(dataset_folder, dataset_settings[args.experiment].mea)

    # init new measure matrix
    if has_prec_mea:
        logger.warning('precomputed measure {} already present.'.format(mea_type))
        if not args.force:
            logger.warning('I quitting the procedure.')
            return
        logger.warning('I am going to recomputed/update it.')
        if mea_type==DISTANCE:
            dataset = DataSet.load_dataset(path=dataset_folder, name=args.experiment,
                                           precomputed_distance=False,
                                           precomputed_kernel=True)
            mat_new = dataset.prec_distance_mat
        elif mea_type==KERNEL:
            dataset = DataSet.load_dataset(path=dataset_folder, name=args.experiment,
                                           precomputed_distance=False,
                                           precomputed_kernel=True)
            mat_new = dataset.prec_kernel_mat
        notes = 'updating precomputed matrix... '
    else:
        notes = 'generated from scratch... '
        dataset = DataSet.load_dataset(path=dataset_folder, name=args.experiment,
                                       precomputed_distance=False,
                                       precomputed_kernel=False)
        mat_new = np.empty((dataset.no_graphs, dataset.no_graphs))
        mat_new[:] = np.nan

    # create measure_matrix
    dataset_settings[args.experiment].mea.set_n_jobs(n_jobs=args.noJobs)
    measure_fun = dataset_settings[args.experiment].mea.get_measure_fun(verbose=True)
    # if dataset_settings[args.experiment].cla is None:
    source_classes = list(dataset.elements.keys())
    target_classes = list(dataset.elements.keys())
    try:
        if dataset_settings[args.experiment].cla is not None:
            raise NotImplementedError
        # source_classes, target_classes = dataset_settings[args.experiment].cla
    except:
        pass
    g_format = dataset_settings[args.experiment].gformat
    source_idx, source = dataset.get_graphs(classes=source_classes, format=['idx', g_format])
    target_idx, target = dataset.get_graphs(classes=target_classes, format=['idx', g_format])

    measure_submat = measure_fun(source, target)

    # update distance matrix
    mat_new[np.array(source_idx)[..., None], np.array(target_idx)[None, ...]] = measure_submat
    today = datetime.now().strftime('%G%m%d%H%M-%f_cpm')
    notes = '{}{}'.format(notes, type(dataset_settings[args.experiment].mea))
    if mea_type==DISTANCE:
        dataset.store(path=dataset_folder, dist_mat=mat_new, dist_gen=today, dist_notes=notes)
    elif mea_type==KERNEL:
        dataset.store(path=dataset_folder, kernel_mat=mat_new, kernel_gen=today, kernel_notes=notes)


if __name__ == "__main__":
    main(sys.argv[1:])
