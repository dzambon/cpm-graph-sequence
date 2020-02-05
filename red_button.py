# Set the debugging options
from cdg.utils import logger, gather_results
logger.set_stdout_level(logger.DEBUG)
rb_logger, _ = logger.create_new_logger(name='red_button',
                                     stdout_level=logger.DEBUG,
                                     filelog_level=logger.DEBUG,
                                     date_in_fname=False)

import argparse
parser = argparse.ArgumentParser(description='')

parser.add_argument('--debug', action='store_true')
parser.add_argument('--debugSBM', action='store_true')

parser.add_argument('--generate', action='store_true')
parser.add_argument('--precompute', action='store_true')

parser.add_argument('--all', action='store_true')
#delaunay
parser.add_argument('--delaunay', action='store_true')
parser.add_argument('--delstd', action='store_true')
parser.add_argument('--delratio', action='store_true')
parser.add_argument('--delmulti', action='store_true')
#iam
parser.add_argument('--iam', action='store_true')
parser.add_argument('--aids', action='store_true')
parser.add_argument('--mutag', action='store_true')
parser.add_argument('--letter', action='store_true')
#kaggle
parser.add_argument('--kaggle', action='store_true')
parser.add_argument('--dog', action='store_true')
parser.add_argument('--human', action='store_true')
#sbm (supplementary material)
parser.add_argument('--supplmat', action='store_true')
parser.add_argument('--sbmnc', action='store_true')
parser.add_argument('--sbmpg', action='store_true')
parser.add_argument('--sbmimf', action='store_true')
#drawplot
parser.add_argument('--figures', action='store_true')
parser.add_argument('--statPvalPlot', action='store_true')
parser.add_argument('--sbmDriftPlot', action='store_true')
parser.add_argument('--delDistPlot', action='store_true')
#results
parser.add_argument('--noresults', action='store_true')
parser.add_argument('--readfrom', type=str, default='./',
                    help='prefix for the folder with the results')

error_list = []

from config import SEED, NO_JOBS, REP, FOLDER_PREFIX, SELECTED_CPM

def main(args):
    rb_logger.info('Red Button has been pressed.')
    
    #Sei lists
    if args.all:
        args.delaunay = True
        args.iam = True
        args.kaggle = True
        args.supplmat = True
    if args.delaunay:
        args.delstd = True
        args.delratio = True
        args.delmulti = True
    if args.iam:
        args.aids = True
        args.mutag = True
        args.letter = True
    if args.kaggle:
        args.dog = True
        args.human = True
    if args.supplmat:
        args.sbmnc = True
        args.sbmpg = True
        args.sbmimf = True
    if args.figures:
        args.statPvalPlot = True
        args.sbmDriftPlot = True
        args.delDistPlot = True

    # Generate Delaunay data set, if missing
    if args.generate:
        import generate_dataset
        rb_logger.info('Generate Delaunay data set.')
        # generate_delaunay_dataset.main(['--complexity'])
        generate_dataset.main(['--complexity'], logger=rb_logger)
        # generate_delaunay_dataset.main(['--variation'])
        generate_dataset.main(['--variation'], logger=rb_logger)
        #generate_dataset.main(['--dcsbmfrag'])

    # Precompute graph measures
    if args.precompute:
        # delaunay
        import precompute_measures
        from settings import dataset_settings
        rb_logger.info('Precompure measures.')
        for d in dataset_settings.keys():
            precompute_measures.main(['--noJobs', str(NO_JOBS), '--experiment', d], logger=rb_logger)

    pars_std = ['--numSimulations', str(REP), '--seed', str(SEED), '--folderPrefix', FOLDER_PREFIX, '--noJobs', str(NO_JOBS)]

    fe = 0  #failed_experiemtns = 0
    # debug
    if args.debug:
        for cpm in get_cpm(['mucpm'], SELECTED_CPM):
            from settings import sbm_stat_list
            for s in sbm_stat_list:
                fe += run_exp(pars_std + ['--test', cpm, '--experiment', s])

    if args.debugSBM:
        for cpm in get_cpm(['mucpm', 'ecpm', 'ediv'], SELECTED_CPM):
            from settings import sbm_debug_list
            for s in sbm_debug_list:
                fe += run_exp(pars_std + ['--test', cpm, '--experiment', s])

    # delaunay
    for cpm in get_cpm(['mucpm', 'ecpm', 'dcpm', 'kcpm', 'ediv'], SELECTED_CPM):
        if args.delstd:
            from settings import delaunay_class1_list
            fe += run_exp(pars_std + ['--test', cpm, '--experiment', 'del', '--particular'] + [str(i) for i in delaunay_class1_list])
        # if args.delno:
        #     run_exp(pars_std + ['--test', cpm, '--experiment', 'del', '--particular', -1])
        if args.delratio:
            from settings import delaunay_ratio_class1_list
            fe += run_exp(pars_std + ['--test', cpm, '--experiment', 'del_ratio', '--particular'] + [str(i) for i in delaunay_ratio_class1_list])
    if args.delmulti and 'ediv' in SELECTED_CPM:
        fe += run_exp(pars_std + ['--test', 'ediv', '--experiment', 'del_multi1'])
        fe += run_exp(pars_std + ['--test', 'ediv', '--experiment', 'del_multi2'])

    # iam
    for cpm in get_cpm(['mucpm', 'ecpm', 'ediv'], SELECTED_CPM):
        if args.aids:
            fe += run_exp(pars_std + ['--test', cpm, '--experiment', 'aids'])
            fe += run_exp(pars_std + ['--test', cpm, '--experiment', 'aids_no'])
        if args.mutag:
            fe += run_exp(pars_std + ['--test', cpm, '--experiment', 'mutag'])
            fe += run_exp(pars_std + ['--test', cpm, '--experiment', 'mutag_no'])
    if args.letter and 'ediv' in SELECTED_CPM:
        fe += run_exp(pars_std + ['--test', 'ediv', '--experiment', 'letter_multi1'])
        fe += run_exp(pars_std + ['--test', 'ediv', '--experiment', 'letter_multi2'])

    # kaggle-seizure
    for cpm in get_cpm(['mucpm', 'ecpm', 'ediv'], SELECTED_CPM):
        if args.dog:
            from settings import kaggle_dog_list
            for dog in kaggle_dog_list:
                fe += run_exp(pars_std + ['--test', cpm, '--experiment', 'dog{}'.format(dog)])
        if args.human:
            from settings import kaggle_human_list
            for hum in kaggle_human_list:
                fe += run_exp(pars_std + ['--test', cpm, '--experiment', 'hum{}'.format(hum)])

    # sbm (suppl mat)
    for cpm in get_cpm(['mucpm', 'ecpm', 'ediv'], SELECTED_CPM):
        from settings import sbm_nc_list, sbm_imf_list, sbm_prop_glob_list
        sbm_to_do = []
        if args.sbmnc:
            sbm_to_do += sbm_nc_list
        if args.sbmimf:
            sbm_to_do += sbm_imf_list
        if args.sbmpg:
            sbm_to_do += sbm_prop_glob_list
        for s in sbm_to_do:
            fe += run_exp(pars_std + ['--test', cpm, '--experiment', s])

    # Plot example
    if args.statPvalPlot:
        for cpm in get_cpm(['mucpm', 'ecpm'], SELECTED_CPM):
            fe += run_exp(['--numSimulations', '1', '--seed', str(SEED),
                           '--folderPrefix', 'tmp_drawplot', '--statPvalPlot',
                           '--test', cpm, '--experiment', 'del', '--particular', '9'])
    if args.sbmDriftPlot:
        from settings import sbm_driftplot_list
        for cpm in get_cpm(['mucpm'], SELECTED_CPM):
            for s in sbm_driftplot_list:
                fe += run_exp(['--numSimulations', '1', '--seed', str(SEED),
                               '--folderPrefix', 'tmp_drawplot', '--sbmDriftPlot',
                               '--test', cpm, '--experiment', s])
    if args.delDistPlot:
        from show_del_distances import main as draw_dist
        draw_dist()

    # Read results
    if not args.noresults:
        from run_simulation import DataSet_dynamicDCSBM, SimulationCPM_dynamicDCSBM        # these import are necessary to deserialize the package
        import os, tempfile
        exp_list = tempfile.NamedTemporaryFile()
        with open(exp_list.name, 'w') as f:
            # This is code is because I'm lazy. I open the file here so that I can simple close it at the end.
            if args.readfrom.endswith('.txt'):
                zipfile_list = '{}/{}'.format(os.getcwd(), args.readfrom)
            else:
                logger.info('generating the list of zipfile...')
                working_dir = '{}/{}'.format(os.getcwd(), args.readfrom)
                for filename in os.listdir(working_dir):
                    if filename.startswith(FOLDER_PREFIX) and filename.endswith('.zip'):
                        # for filename in os.listdir(FOLDER_PREFIX+'*.zip'):
                        #exp_list.write(bytes('{}/{}\n'.format(working_dir, filename), 'utf-8'))
                        f.writelines('{}/{}\n'.format(working_dir, filename))
                        #print('{}/{}'.format(working_dir, filename))
                        exp_list.seek(0)
                        f.flush()
                # find $PWD/$FOLDER_PREFIX*.zip > exp_list.txt
                #print(exp_list.name)
                zipfile_list = exp_list.name
            table = gather_results(args_to_parse=['-f', zipfile_list, '-l', 'tsp'], setting_list=['dataset', 'classes', 'cpm'])
            f.close()
            exp_list.close()
        
        print(post_process_table(table))
        
    if fe > 0:
        import time
        time.sleep(3)
        rb_logger.warning("Number of failed experiments: {}".format(fe))

# Run actual experiments
def run_exp(arg_list):
    import run_simulation
    global rb_logger, error_list
    try:
        rb_logger.info('running {}'.format(arg_list))
        run_simulation.main(arg_list)
        rb_logger.info('concluded.')
        return 0
    except Exception:
        import sys, traceback
        exc_type, exc_value, exc_tb = sys.exc_info()
        exc_str = ' '.join(traceback.format_exception(exc_type, exc_value, exc_tb))
        traceback.print_exception(exc_type, exc_value, exc_tb)
        rb_logger.error(exc_str)
        error_list.append(exc_str)
        return 1


def post_process_table(latex_table):
    import re
    # See also: 
    #  - gather results in cdg.utils.extra 
    #  - https://www.w3schools.com/python/python_regex.asp

    line_header = 3
    latex_split = latex_table.split('\n')
    
    #header
    h_split = ['experiment'] + latex_split[line_header].split('&')[3:]
    h_split[-1] = h_split[-1][:-2]
    h_name = ['\\textbf{{{}}}'.format(h.replace(' ', '')) for h in h_split]
    header = ''
    for h in h_name:
        header += ' {} &'.format(h)
    latex_split[line_header] = header[:-1] + '\\\\'
    
    for i in range(line_header+1, len(latex_split)):
        # Parse method
        latex_split[i] = re.sub(r'CPM\+t-test\s*', r'\\mucpm ', re.escape(latex_split[i]))
        latex_split[i] = re.sub(r'ECPM\s*', r'\\ecpm  ', latex_split[i])
        latex_split[i] = re.sub(r'EDivR\s*', r'\\ediv  ', latex_split[i])
        latex_split[i] = re.sub(r'EDivR\s*', r'\\ediv  ', latex_split[i])
        latex_split[i] = re.sub(r'mmdCPM\s*', r'\\kcpm  ', latex_split[i])
        latex_split[i] = re.sub(r'DCPM\s*', r'\\dcpm  ', latex_split[i])

        # Parse Delaunay
        latex_split[i] = re.sub(r'd\.Del\s*\&\s*Delaunay\s*\&\s((\.\d*)*)\s*\&',      r'Del\\textsubscript{\1}       &', latex_split[i])
        latex_split[i] = re.sub(r'd\.Del\s*\&\s*DelaunayMulti\s*\&\s((\.\d*)*)\s*\&', r'Del\\textsubscript{\1}       &', latex_split[i])
        latex_split[i] = re.sub(r'd\.Del\s*\&\s*DelaunayRatio\s*\&\s((\.\d*)*)\s*\&', r'Del*\\textsubscript{\1}      &', latex_split[i])
        # Parse Mutag classes
        latex_split[i] = re.sub(r'\.nonmutagen([\.,\s])', r'.0\1', latex_split[i])
        latex_split[i] = re.sub(r'\.mutagen([\.,\s])', r'.1\1', latex_split[i])
        # Parse AIDS classes
        latex_split[i] = re.sub(r'\.i([\.,\s])', r'.0\1', latex_split[i])
        latex_split[i] = re.sub(r'\.a([\.,\s])', r'.1\1', latex_split[i])
        # Parse Kaggle classes
        latex_split[i] = re.sub(r'\.preictal([\.,\s])', r'.0\1', latex_split[i])
        latex_split[i] = re.sub(r'\.ictal([\.,\s])', r'.1\1', latex_split[i])
        # Parse IAM id 
        latex_split[i] = re.sub(r'd\.Mut\s*\&[\s,a-z,A-Z]*\&\s((\.\d)*)\s*\&',    r'Mut\\textsubscript{\1}       &', latex_split[i])
        latex_split[i] = re.sub(r'd\.AIDS\s*\&[\s,a-z,A-Z]*\&\s((\.\d)*)\s*\&',   r'AIDS\\textsubscript{\1}      &', latex_split[i])
        latex_split[i] = re.sub(r'd\.Let\s*\&[\s,a-z,A-Z]*\&\s((\.[A-Z])*)\s*\&', r'Let\\textsubscript{\1}       &', latex_split[i])
        # Parse Kaggle id
        latex_split[i] = re.sub(r'd\.Hum\s*\&\s*Human(\d*)\s*\&\s((\.\d)*)\s*\&', r'H\1\\textsubscript{\2}       &', latex_split[i])
        latex_split[i] = re.sub(r'd\.Dog\s*\&\s*Dog(\d*)\s*\&\s((\.\d)*)\s*\&',   r'D\1\\textsubscript{\2}       &', latex_split[i])
        # Clean class name
        latex_split[i] = re.sub(r'textsubscript{\.([a-z,A-Z,\d])', r'textsubscript{\1', latex_split[i])

        # Parse SBM
        latex_split[i] = re.sub(r'd\.SBM\s*\&\s*SBMsp-globa\.(\d*)\%v\.(\d*)\s*\&\s*\.None\s*\&', r'Global  & \1 & \2 &', latex_split[i])
        latex_split[i] = re.sub(r'd\.SBM\s*\&\s*SBMsp-propa\.(\d*)\%v\.(\d*)\s*\&\s*\.None\s*\&',  r'Prop    & \1 & \2 &', latex_split[i])
        latex_split[i] = re.sub(r'd\.SBM\s*\&\s*SBM-i1090a\.(\d*)\%v\.(\d*)\s*\&\s*\.None\s*\&',   r'Intense & \1 & \2 &', latex_split[i])
        latex_split[i] = re.sub(r'd\.SBM\s*\&\s*SBM-m5050a\.(\d*)\%v\.(\d*)\s*\&\s*\.None\s*\&',   r'Merge   & \1 & \2 &', latex_split[i])
        latex_split[i] = re.sub(r'd\.SBM\s*\&\s*SBM-f5050a\.(\d*)\%v\.(\d*)\s*\&\s*\.None\s*\&',   r'Frag    & \1 & \2 &', latex_split[i])

        while True:
            tmp = re.sub(r'subscript{([a-z,A-Z,\d,\,]*)\.([a-z,A-Z,\d,\,]*)', r'subscript{\1,\2', latex_split[i])
            if tmp == latex_split[i]:
                break
            latex_split[i] = tmp

    latex_table_new = ''
    for l in latex_split:
        latex_table_new += l + "\n"

    return latex_table_new

def get_cpm(list1, list2):
    if list2 is None:
        return list(set(list1))
    else: 
        return list(set(list1) & set(list2)) 
    
if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
