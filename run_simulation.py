# Set the debugging options
from cdg.utils import logger
try:
    logger.set_stdout_level(logger.DEBUG)
    log_file, idcode = logger.set_filelog_level(level=logger.DEBUG)
    logger.enable_logrun(level=False)
except:
    pass

import sys
import argparse
import numpy as np

from cdg import SimulationCPM, DataSet
from cdg.simulation import ParametersCPM, create_foldername

from config import SIGNIFICANCE_LEVEL
from settings import CPM_TEST_INSTANCES, LIST_CPM, MUCPM, ECPM, DCPM, KCPM, EDIV
from settings import experiment_settings, experiments_with_multiple_settings
from settings import EMBEDDINGS, NO_EMBEDDING_USE_DISTANCE, NO_EMBEDDING_USE_KERNEL

# List of all possible experiments
list_all_experiments = [e.id for e in experiment_settings]
# List of available experiments per CPM test
_available_mucpm = {}
for e in list_all_experiments:
    _available_mucpm[e] = False
_available_ecpm = _available_mucpm.copy()
_available_dcpm = _available_mucpm.copy()
_available_kcpm = _available_mucpm.copy()
_available_ediv = _available_mucpm.copy()
# check that are actually possible
for k in _available_mucpm: assert k in list_all_experiments
for k in _available_ecpm:  assert k in list_all_experiments
for k in _available_dcpm:  assert k in list_all_experiments
for k in _available_kcpm:  assert k in list_all_experiments
for k in _available_ediv:  assert k in list_all_experiments

# Define the argument parser
parser = argparse.ArgumentParser(description='')
# parser.add_argument('-e', '--experiment', type=str, help='experiment name')
parser.add_argument('-e', '--experiment', type=str,
                    required=True,
                    choices=list_all_experiments,
                    help='experiment name')
parser.add_argument('-p', '--particular', nargs='+', type=int,
                    default=True,
                    help='particular type of that experiment')
parser.add_argument('-t', '--test', type=str,
                    required=True,
                    choices=[MUCPM, ECPM, EDIV, DCPM, KCPM],
                    help='tyep of test')
parser.add_argument('-r', '--numSimulations', type=int, default=15,
                    help='number of repreated simulations')
parser.add_argument('-s', '--seed', type=int, help='seed for the PRG')
parser.add_argument('-f', '--folderPrefix', type=str, default='cpm1_',
                    help='prefix for the folder with the results')
parser.add_argument('-j', '--noJobs', type=int, default=1,
                    help='jobs in joblib')
parser.add_argument('--statPvalPlot', action='store_true',
                    help='whether to plot the cpm stats and pvals')
parser.add_argument('--sbmDriftPlot', action='store_true',
                    help='whether to plot the drift in the DCSBM sequences')


bluecustom = (.05, .4, .6)
def draw_statPvalPlot(stats, pvals, cp_tr, cp_est, name):
    """
    Illustrative example of how the CPM works. 
    It draws the statistics of each two-sample test and the associated p-values. 
    """
    cp_est = np.nanargmax(stats)
    margin = 15
    xticks = [i * 10 for i in range(len(stats)//10)] + [len(stats)]
    xticks = [margin, len(stats)-margin]

    pvals[np.where(pvals == 2.)] = np.nan

    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import host_subplot
    import matplotlib.ticker as ticker
    from matplotlib import rc
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    # import mpl_toolkits.axisartist as AA

    plt.figure(figsize=(5, 2))
    ax1 = host_subplot(111)
    # fig, ax1 = plt.subplots()
    plt.subplots_adjust(right=.92)
    plt.subplots_adjust(left=0.1)
    plt.subplots_adjust(bottom=0.20)
    plt.subplots_adjust(top=.85)
    ax1.set_xlabel("time index")
    ax1.set_ylabel("$p$-value")

    ax2 = ax1.twinx()
    ax2.set_ylabel("Statistic $s_e(t)$", color=bluecustom)

    ax1.set_xticks(np.array(xticks + [cp_tr, cp_est]))
    #ax1.set_xticklabels([str(l) for l in xticks] + ['$t*$', '$\hat t$'])

    # par1.plot(time_steps, pval_g, 'r:', label="graph p-value")

    if cp_tr == cp_est:
        print('tik==')
        ax1.set_xticklabels(['margin', 'T-margin', '$t*=\hat t$'])
    else:
        print('tik!=')
        ax1.set_xticklabels(['margin', 'T-margin', '$t*$', '$\hat t$'])

    #ax1.set_xticklabels(['margin', 'T-margin', '$t*$', '$\hat t$'])
    time_steps = [i for i in range(len(stats))]
    ax2.plot(time_steps, stats, '--', label="$s_e(t)$", color=bluecustom)
    # host.plot(time_steps, th, 'k--', label="threshold")
    ax1.semilogy(time_steps, pvals, 'k', label="p-value")

    ytick_pval = [0.01, 0.05, 1]
    #ax2.set_yscale('log')
    ax1.get_yaxis().set_major_formatter(ticker.ScalarFormatter())
    ax1.set_yticks(ytick_pval)
    ax2.tick_params('y', colors=bluecustom)
    ax1.grid(True)

    ax1.legend(loc='lower right')
    plt.draw()
    plt.savefig('cpm_example_{}.pdf'.format(name))
    plt.show()


class DataSet_dynamicDCSBM(DataSet):
    """ 
    Custom sub-class that deals with the generation of the DCSBM graphs on the fly. 
    """

    def __init__(self, setting=None, no_nodes=None, alpha_update=None, **kwargs):
        self.setting = setting
        self.no_nodes = no_nodes
        self.alpha_update = alpha_update
        super().__init__(**kwargs)


class SimulationCPM_distance(SimulationCPM):
    """ 
    Custom sub-class that bypass the embedding stage and perform a distance CPM
    directly on the sequence of graphs. 
    """
    
    def run_single_simulation(self, **kwargs):
        self.log.debug('sequence generation...')
        return_indices = self.dataset.has_prec_distance()
        g_train, g_test, change_points = self.sequence_generator(dataset=self.dataset, pars=self.pars,
                                                                 return_indices=return_indices)
        self.changes_true.append(change_points)
        
        self.log.debug('(dis)similarity matrix computation...')
        dm = self.dataset.distance_measure(g_test, g_test)
        
        self.log.debug('operating phase...')
        self.cpm.reset()
        assert self.use_fwer
        change_points_est, _ = self.cpm.predict(dm, is_dist_mat=True, alpha=self.pars.significance_level,
                                                margin=self.pars.margin, fwer=self.use_fwer, **kwargs)
        self.changes_est.append(change_points_est)
        
        return change_points_est


class SimulationCPM_kernel(SimulationCPM):
    """ 
    Custom sub-class that bypass the embedding stage and perform a kernel CPM
    directly on the sequence of graphs. 
    """
    
    def run_single_simulation(self, **kwargs):
        self.log.debug('sequence generation...')
        return_indices = self.dataset.has_prec_kernel()
        g_train, g_test, change_points = self.sequence_generator(dataset=self.dataset, pars=self.pars,
                                                                 return_indices=return_indices)
        self.changes_true.append(change_points)
        
        self.log.debug('(dis)similarity matrix computation...')
        km = self.dataset.kernel_measure(g_test, g_test)
        
        self.log.debug('operating phase...')
        self.cpm.reset()
        assert self.use_fwer
        change_points_est, _ = self.cpm.predict(km, is_kernel_mat=True, alpha=self.pars.significance_level,
                                                margin=self.pars.margin, fwer=self.use_fwer, **kwargs)
        self.changes_est.append(change_points_est)
        
        return change_points_est


class SimulationCPM_dynamicDCSBM(SimulationCPM):
    """ Custom sub-class that generates the considered DCSBM on the fly. """
    
    def __init__(self, cpm=None, setting=None, no_nodes=None):
        super().__init__(cpm=cpm)
        self.setting = setting
        self.no_nodes = no_nodes
    
    @staticmethod
    def get_DCSBM_paramters(setting, no_nodes):
        # Figure 2 wilson et al
        
        if setting[:4] == 'glob':
            comm0, prob_mat0, delta0 = [list(range(no_nodes))], [[.2]], [.5]
            comm1, prob_mat1, delta1 = comm0, [[.25]], delta0
            distrib = 'poisson'
        
        elif setting[:5] == 'i1090':
            no_el_c0 = int(no_nodes * .1)
            comm0 = [list(range(no_el_c0)), list(range(no_el_c0, no_nodes))]
            prob_mat0 = [[.3, .1], [.1, .3]]
            delta0 = [.5] * 2
            comm1, prob_mat1, delta1 = comm0, [[.4, .1], [.1, .3]], delta0
            distrib = 'poisson'
        
        elif setting[:5] == 'm5050':
            no_el_c0 = no_nodes // 2
            comm0 = [list(range(no_el_c0)), list(range(no_el_c0, no_nodes))]
            prob_mat0 = [[.3, .1], [.1, .3]]
            delta0 = [.5] * 2
            comm1, prob_mat1, delta1 = comm0, [[.2, .2], [.2, .2]], delta0
            distrib = 'poisson'
        
        elif setting[:5] == 'f5050':
            no_el_c0 = no_nodes // 2
            comm0 = [list(range(no_el_c0)), list(range(no_el_c0, no_nodes))]
            prob_mat0 = [[.3, .1], [.1, .3]]
            delta0 = [.5] * 2
            comm1, prob_mat1, delta1 = comm0, [[.4, .2], [.2, .1]], delta0
            distrib = 'poisson'
        
        elif setting[:5] == 'prop':
            comm0, prob_mat0, delta0 = [list(range(no_nodes))], [[.2]], [(.5, .5)]
            comm1, prob_mat1, delta1 = comm0, [[.2]], [(.5, 2.)]
            distrib = 'poisson'
        
        elif setting[:8] == 'noch5050':
            # No change sequence. I had to set delta to zero, so that I will have
            # the same theta. This was necessary to keep the same framework in
            # which I generate two sequences from two models.
            no_el_c0 = no_nodes // 2
            comm0 = [list(range(no_el_c0)), list(range(no_el_c0, no_nodes))]
            prob_mat0 = [[.3, .1], [.1, .3]]
            delta0 = [.0] * 2
            comm1, prob_mat1, delta1 = comm0, prob_mat0, delta0
            distrib = 'poisson'

        elif setting[:8] == 'gllowden':
            # global change about 20% density and delta=0
            comm0, prob_mat0, delta0 = [list(range(no_nodes))], [[.2]], [0.]
            comm1, prob_mat1, delta1 = comm0, [[.25]], delta0
            distrib = 'uniform'
        
        elif setting[:8] == 'glmidden':
            # global change about 50% density and delta=0
            comm0, prob_mat0, delta0 = [list(range(no_nodes))], [[.5]], [0.]
            comm1, prob_mat1, delta1 = comm0, [[.55]], delta0
            distrib = 'uniform'
        
        elif setting[:9] == 'glhighden':
            # global change about 80% density and delta=0
            comm0, prob_mat0, delta0 = [list(range(no_nodes))], [[.8]], [0.]
            comm1, prob_mat1, delta1 = comm0, [[.75]], delta0
            distrib = 'uniform'


        elif setting[:12] == 'gllowden_poi':
            # global change about 20% density and delta=0
            comm0, prob_mat0, delta0 = [list(range(no_nodes))], [[.2]], [0.]
            comm1, prob_mat1, delta1 = comm0, [[.25]], delta0
            distrib = 'poisson'

        elif setting[:12] == 'glmidden_poi':
            # global change about 50% density and delta=0
            comm0, prob_mat0, delta0 = [list(range(no_nodes))], [[.5]], [0.]
            comm1, prob_mat1, delta1 = comm0, [[.55]], delta0
            distrib = 'poisson'

        elif setting[:13] == 'glhighden_poi':
            # global change about 80% density and delta=0
            comm0, prob_mat0, delta0 = [list(range(no_nodes))], [[.8]], [0.]
            comm1, prob_mat1, delta1 = comm0, [[.75]], delta0
            distrib = 'poisson'
            
        else:
            raise ValueError('setting {} not recognized'.format(setting))
        
        return comm0, prob_mat0, delta0, comm1, prob_mat1, delta1, distrib
    
    @classmethod
    def sequence_generator(cls, dataset, pars, return_indices=True):
        # parse possible None parameters
        assert len(pars.subseq_lengths_t) == 2
        
        from cdg.graph import DegreeCorrectedStochasticBlockModel, DynamicsGenerator
        
        comm0, prob_mat0, delta0, comm1, prob_mat1, delta1, distrib = cls.get_DCSBM_paramters(dataset.setting,
                                                                                              dataset.no_nodes)
        model0 = DegreeCorrectedStochasticBlockModel(communities=comm0, prob_matrix=prob_mat0, delta=delta0)
        model1 = DegreeCorrectedStochasticBlockModel(communities=comm1, prob_matrix=prob_mat1, delta=delta1)

        dyn_model0 = DynamicsGenerator(alpha=dataset.alpha_update,
                                       getter=lambda: model0.get(distrib=distrib, no_graphs=1)[0])
        dyn_model1 = DynamicsGenerator(alpha=dataset.alpha_update,
                                       getter=lambda: model1.get(distrib=distrib, no_graphs=1)[0])
        
        G0 = dyn_model0.get(no_graphs=pars.subseq_lengths_t[0])
        G1 = dyn_model1.get(no_graphs=pars.subseq_lengths_t[1] + 1, graph_seed=G0[-1])[1:]
        
        tr_len0 = int(pars.train_len_ratio * pars.subseq_lengths_t[0])
        tr_len1 = int(pars.train_len_ratio * pars.subseq_lengths_t[1])
        
        # split train an test set randombly
        p0 = np.random.permutation(pars.subseq_lengths_t[0])
        p1 = np.random.permutation(pars.subseq_lengths_t[1])
        p0tr = sorted(p0[:tr_len0])
        p1tr = sorted(p1[:tr_len1])
        p0te = sorted(p0[tr_len0:])
        p1te = sorted(p1[tr_len1:])
        
        g_train = [G0[i] for i in p0tr] + [G1[i] for i in p1tr]
        g_test  = [G0[i] for i in p0te] + [G1[i] for i in p1te]

        if pars.subseq_lengths_t[0] == 0 or pars.subseq_lengths_t[1] == 0:
            changes = []
        else:
            changes = [pars.subseq_lengths_t[0] - tr_len0]

        return g_train, g_test, changes

class SimulationCPM_dynamicDCSBM_drifts(SimulationCPM_dynamicDCSBM):
    """ Custom sub-class to create figures with the embedding vectors. """

    def run_single_simulation(self, pass_to_next=None, **kwargs):
    
        self.log.debug('sequence generation...')
        return_indices = self.dataset.has_prec_distance()
        g_train, g_test, change_points = self.sequence_generator(dataset=self.dataset, pars=self.pars,
                                                                 return_indices=return_indices)
        self.changes_true.append(change_points)
    
        self.log.debug('embedding...')
        self.pars.embedding_method.fit(graphs=g_train, dist_fun=self.dataset.distance_measure, **kwargs)
        x = self.pars.embedding_method.transform(data=g_test)
    
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib import rc
        # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        # rc('font',**{'family':'serif','serif':['Palatino']})
        rc('text', usetex=True)
        T = x.shape[0]
        c = cm.rainbow(np.linspace(0, 1, T))

        fig = plt.figure(figsize=(5,3.5))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x[:, 0], x[:, 1], x[:, 2], color=c)
        ax.set_xlabel(r"$x_0$", fontsize=15)
        ax.set_ylabel(r"$x_1$", fontsize=15)
        ax.set_zlabel(r"$x_2$", fontsize=15)
        plt.title(r"$\rho$ = {:0.1f}".format(self.dataset.alpha_update), fontsize=18)
        # plt.tight_layout(pad=0)
        plt.draw()
        plt.savefig("nc_{}.pdf".format(self.dataset.name.replace("%", "_")))
        # plt.show()

        fig = plt.figure(figsize=(5, 1.4))
        ax = fig.add_subplot(111)
        ax.scatter(np.arange(T), np.zeros(T), color=c)
        ax.set_xlabel(r"time index t")
        ax.set_yticks([])
        ax.set_xticks([0, T-1])
        ax.set_xticklabels(["0", "T"])
        plt.title(r"Color map", fontsize=15)
        plt.subplots_adjust(left=0.03, right=0.97, top=0.8, bottom=0.35)
        plt.draw()
        plt.savefig("nc_colormap.pdf".format(self.dataset.name.replace("%", "_")))

        self.changes_est.append([])
        return []


class SimulationCPM_with_plot(SimulationCPM):
    """ Custom sub-class to create figures with obtained statistics and p-values. """
    def run_single_simulation(self, pass_to_next=None, **kwargs):
        changes = super().run_single_simulation(pass_to_next=pass_to_next)
        if len(changes)==0:
            changes.append(-1)
        draw_statPvalPlot(stats=self.cpm.stats_seq, pvals=self.cpm.pvals_fwer_seq,
                          cp_tr=self.changes_true[0][0], cp_est=self.cpm.cps_fwer[0],
                          name=self.cpm.name)
        self.log.info('{}) cps {}\t p_fwer {}'.format(self.cpm.name, self.cpm.cps_fwer, self.cpm.pvals_fwer))
        return changes.pop(-1)


def run(cpm, pars, path, name, repetitions, seed, prefix, dyn_setup, n_jobs, plot_res=False, plot_drift=False):
    '''
    aux function that run the actual experiment
    '''
    #create output folder
    folder = create_foldername(prefix=prefix, exp_name=name, pars=pars, cpm=cpm) + idcode

    #identify the correct simulation
    use_distance = True
    if plot_res:
        cpm._use_max_stat = False
        simulation = SimulationCPM_with_plot(cpm=cpm)
    elif plot_drift:
        simulation = SimulationCPM_dynamicDCSBM_drifts(cpm=cpm)
    elif dyn_setup is not None:
        if pars.embedding_method == NO_EMBEDDING_USE_DISTANCE(0) \
                or pars.embedding_method == NO_EMBEDDING_USE_KERNEL(0):
            raise NotImplementedError()
        simulation = SimulationCPM_dynamicDCSBM(cpm=cpm)
    elif pars.embedding_method == NO_EMBEDDING_USE_DISTANCE(0):
        simulation = SimulationCPM_distance(cpm=cpm)
    elif pars.embedding_method == NO_EMBEDDING_USE_KERNEL(0):
        use_distance = False
        simulation = SimulationCPM_kernel(cpm=cpm)
    else:
        simulation = SimulationCPM(cpm=cpm)

    #load dataset
    if dyn_setup is not None:
        dataset = DataSet_dynamicDCSBM(name=name, setting=dyn_setup['setting'],
                                       no_nodes=dyn_setup['no_nodes'], alpha_update=dyn_setup['alpha_update'],
                                       distance_measure=dyn_setup['distance'])
        pars.subseq_lengths_t = dyn_setup['seq_lengths']
    else:
        dataset = DataSet.load_dataset(path=path, name=name,
                                       precomputed_distance=use_distance, precomputed_kernel=not use_distance,
                                       skip_graphs=True)
        pars.subseq_lengths_t = [len(dataset.elements[c]) for c in pars.classes]
    # if pars.train_len_t > 300:
    #     pars.train_len_t = 300
    pars.freeze()

    #set up simulation
    simulation.set(parameters=pars, dataset=dataset, no_simulations=repetitions, folder='./'+folder)

    #run the simulation
    logger.info('    ****     ****    ****')
    logger.info('running {} with {}'.format(name, cpm.name))
    simulation.run(seed=seed, logfile=log_file, n_jobs=n_jobs)
    logger.info('completed run of {} with {}'.format(dataset.name, cpm.name))
    # os.remove(log_file)
    # open(log_file, 'a').close()


def get_single_exp_setting(exp_list):
    '''
    Select the settings of the experiment to perform.
    The idea is to prune all the unwanted experiments.
    :param exp_list:
    :param dataset_path:
    :return:
    '''
    
    # init outputs
    pars, path, name, dyn_setup = None, None, None, None
    exp_found = False
    
    # browse all possible experiments and prune the unwanted ones.
    for expset in experiment_settings:
        el = exp_list.get(expset.id, False)
        if el == True or isinstance(el, list):
            classes = expset.cls.copy()
            path = expset.folder
            name = expset.name
            if not isinstance(el, list):
                exp_list.pop(expset.id)
            else:
                # c, exp_list[expset[EXP_ID]] = get_list_el(el)
                c = exp_list[expset.id].pop()
                if c != -1:
                    classes.append(c)
                if len(exp_list[expset.id]) == 0:
                    exp_list[expset.id] = False
            pars = ParametersCPM()
            pars.classes = classes
            pars.subseq_ratios = [expset.ratio for _ in classes]
            exp_found = True
            dyn_setup = expset.dynamic_setup
            break
        else:
            # prune
            exp_list.pop(expset.id, False)
    
    return path, name, pars, dyn_setup, exp_found


def main(argv):
    """ Identifies the experiment and call the run procedure. """
    logger.info('New call to main: {}'.format(str(argv)))
    args = parser.parse_args(argv)
    # seed
    logger.info('seed set to {}'.format(args.seed))
    # this trick allows me to remove the
    enabled_exp = {MUCPM: _available_mucpm.copy(),
                   ECPM: _available_ecpm.copy(),
                   DCPM: _available_dcpm.copy(),
                   KCPM: _available_kcpm.copy(),
                   EDIV: _available_ediv.copy()}
    
    # Select experiment
    assert args.test in LIST_CPM
    # set_exp(enabled_exp[args.test], args.experiment, args.particular)
    assert args.experiment in enabled_exp[args.test], 'experiment \'{}\' not available'.format(args.experiment)
    enabled_exp[args.test][args.experiment] = args.particular
    if args.experiment in experiments_with_multiple_settings:
        assert args.experiment != True, 'experiment \'{}\' requires a parameter -p (--particular)'.format(args.test)
    
    # browse all tests
    for test_i in LIST_CPM:
        # browse all experiments associated with the test
        while len(enabled_exp[test_i]) > 0:
            # remove from enabled_exp all non disabled experiments, and return the first available
            path, name, pars, dyn_setup, exp_found = get_single_exp_setting(enabled_exp[test_i])
            
            # Perform the experiment
            if exp_found:
                pars.significance_level = SIGNIFICANCE_LEVEL
                pars.embedding_method = EMBEDDINGS[test_i](emb_dim=pars.embedding_dimension)
                run(cpm=CPM_TEST_INSTANCES[test_i], pars=pars, path=path, name=name, repetitions=args.numSimulations,
                    seed=args.seed, prefix=args.folderPrefix, plot_res=args.statPvalPlot, plot_drift=args.sbmDriftPlot,
                    dyn_setup=dyn_setup, n_jobs=args.noJobs)
            
            logger.info('\n\n*******cpm {} : {}'.format(test_i, str(enabled_exp[test_i])))


if __name__ == "__main__":
    main(sys.argv[1:])
