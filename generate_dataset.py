import sys
import argparse
import numpy as np
from collections import OrderedDict

from cdg.utils import logger as cdg_logger
from cdg.graph import DelaunayGraphs, DataSet, isdataset, StochasticBlockModel

from settings import  delaunay_folder, delaunay_var_folder

def delaunay_complexity(args, logger):

    # Parameters
    no_nodes = 7
    no_graphs_per_class = 100
    classes = [0, 8, 9, 10, 11, 12]
    path = delaunay_folder
    
    # Identify if already present
    if isdataset(path):
        logger.warning('Delaunay MuVar data set at path {} is already present'.format(path))
        if args.force:
            logger.warning('Regenerating it...')
        else:
            logger.warning('Quitting the procedure...')
            return
        
    # Create dataset
    graph_dict = DelaunayGraphs().get(seed_points=no_nodes, classes=classes, no_graphs=no_graphs_per_class, sigma=.3)
    gg, y = [], []
    for c in classes:
        gg += graph_dict[c]
        y += [c]*len(graph_dict[c])
        
    return DataSet(graphs=gg, labels=y, name='DelaunayComlexity').store(path=path, notes='label number indicate the complexity (~shift in the mean proportional to (2/3)^c).')

def delaunay_variation(args, logger):

    # Parameters
    no_nodes = 7
    no_graphs_per_class = 100 
    path = delaunay_var_folder

    # Identify if already present
    if isdataset(path):
        logger.warning('Delaunay MuVar data set at path {} is already present'.format(path))
        if args.force:
            logger.warning('Regenerating it...')
        else:
            logger.warning('Quitting the procedure...')
            return
    
    # Create dataset
    sigma0 = .3
    sigmas = [.5, .4, .2, .15]
    del_gen = DelaunayGraphs()
    graph_dict = OrderedDict()
    gg = del_gen.get(classes=[0], no_graphs=no_graphs_per_class, sigma=sigma0,
                     include_seed_graph=True, seed_points=no_nodes)[0]
    graph_dict[int(np.round(sigma0*100))] = gg

    for s in sigmas:
        gg_var = del_gen.get(classes=[0], no_graphs=no_graphs_per_class, sigma=s,
                             include_seed_graph=True, seed_points=del_gen.seed_points)[0]
        graph_dict[int(np.round(s*100))] = gg_var

    gg, y = [], []
    for c in graph_dict.keys():
        print(c)
        gg += graph_dict[c]
        y += [c]*len(graph_dict[c])

    return DataSet(graphs=gg, labels=y, name='DelaunayVar').store(path=path, notes='label number indicate 100*sigma.')

# Define the argument parser
parser = argparse.ArgumentParser(description='')
parser.add_argument('--force', action='store_true', default=False,
                    help='if dataset is already present, this will regenerate from scratch.')
parser.add_argument('--complexity', action='store_true', default=False,
                    help='generate Delaunay data set with increasing complexities.')
parser.add_argument('--variation', action='store_true', default=False,
                    help='generate Delaunay data set with increasing Frèchet variation.')

def main(argv, logger=None):
    if logger is None:
        logger = cdg_logger
        logger.set_stdout_level(logger.DEBUG)
        logger.set_filelog_level(logger.DEBUG)
    args = parser.parse_args(argv)
    if args.complexity:
        delaunay_complexity(args, logger)
    if args.variation:
        delaunay_variation(args, logger)

if __name__ == '__main__':
    main(sys.argv[1:])
