import atlas as basicatlas
from config import settings
import griddef
import glob, os, shutil
import numpy as np
import tqdm

import pick

# Main menu options
menu = [
    ['Clean runs and scripts', 'clean'],
    ['Clean results database', 'clean_results'],
    ['Opacity Distribution Functions', 'dfsynthe'],
    ['Structures', 'atlas'],
    ['Delete failed structures', 'delete_atlas'],
    ['Reindex structures', 'reindex_atlas'],
    ['Line list', 'linelist'],
    ['Spectra', 'rescalc'],
    ['Prepare calculation', 'run'],
]

# Create directories for scripts, runs and results
if not os.path.isdir(settings['scripts_dir']):
    os.mkdir(settings['scripts_dir'])
if not os.path.isdir(settings['runs_dir']):
    os.mkdir(settings['runs_dir'])
if not os.path.isdir(settings['results_dir']):
    os.mkdir(settings['results_dir'])

# Anything to clean?
if len(glob.glob('{}/*'.format(settings['scripts_dir'])) + glob.glob('{}/*'.format(settings['runs_dir']))) > 0:
    menu[0][0] = '** ' + menu[0][0] + ' **'
if len(glob.glob('{}/*'.format(settings['results_dir']))) > 0:
    menu[1][0] = '** ' + menu[1][0] + ' **'

# What are we doing?
action = int(pick.pick([menu[i][0] for i in range(len(menu))], indicator = '=>')[1])


#############################################################################################


dfsynthe_template = """
import atlas

settings = atlas.Settings()
settings.Y = {Y}
settings.zscale = {zscale}
settings.abun = {abun}
atlas.dfsynthe('{output}', settings = settings)
""".strip()

atlas_template = """
import os
import numpy as np
import atlas
import shutil

settings = atlas.Settings()
settings.Y = {Y}
settings.zscale = {zscale}
settings.abun = {abun}
settings.logg = {logg}
restarts = {restarts}
output = {output}
for teff in {teff}:
    settings.teff = teff
    model_name = output[teff]
    restart = restarts[teff]

    # Do not restart from poorly converged models
    if restarts[teff] != 'auto':
        if os.path.isfile(restart + '/failed'):
            restart = 'auto'
        elif np.max(np.abs(atlas.read_structure(restart)[0]['flux_error'])) > 1000:
            restart = 'auto'

    print(' >>> Calculating {{}} from {{}} <<< '.format(model_name, restart))
    try:
        atlas.atlas(model_name, settings = settings, ODF = '{ODF}', restart = restart)
    except Exception as e:
        if e.args[0].find('diverged immediately') != -1:
            shutil.rmtree(model_name)
            print('>>> Model diverged immediately, will try restarting from auto... <<<')
            atlas.atlas(model_name, settings = settings, ODF = '{ODF}', restart = 'auto')
        else:
            raise
    os.remove('{{}}/odf_1.ros'.format(output[teff]))
    os.remove('{{}}/odf_9.bdf'.format(output[teff]))

    # Test the model for integrity and determine convergence
    message = model_name
    while message[-1] == '/':
        message = message[:-1]
    message = os.path.basename(message)
    try:
        structure = atlas.read_structure(model_name)[0]
        e1 = np.max(np.abs(structure['flux_error']))
        e2 = np.max(np.abs(structure['flux_error_derivative']))
        message += '_{{}}_{{}}'.format(e1, e2)
    except:
        message += '_failed'
    f = open('{results_dir}/{{}}'.format(message), 'w')
    f.close()
""".strip()

reindex_atlas_template = """
import os
import atlas
import numpy as np

results_dir = '{results_dir}'
ATLAS_dir = '{ATLAS_dir}'
prefix = '{prefix}'

# Clear existing results
results = {results}
models = {ATLAS}
for result in results:
    os.remove(results_dir + '/' + result)

# Reindexing
for model in models:
    message = model
    model = ATLAS_dir + '/' + model
    if os.path.isfile(model + '/failed'):
        message += '_failed'
    else:
        try:
            meta = atlas.meta(model)
            structure = atlas.read_structure(model)[0]
            e1 = np.max(np.abs(structure['flux_error']))
            e2 = np.max(np.abs(structure['flux_error_derivative']))
            message += '_{{}}_{{}}'.format(e1, e2)
        except:
            message += '_failed'
    f = open('{{}}/{{}}'.format(results_dir, message), 'w')
    f.close()

print('Done {{}}'.format(prefix), flush = True)
""".strip()

rescalc_template = """
import numpy as np
import os, shutil, glob
import atlas
import copy

import sys
sys.path.append('{rescalc_path}')
import rescalc
settings = {settings}
output_dir = '{output_dir}'
results_dir = '{results_dir}'


for output in {output}:

    rescalc.settings = copy.deepcopy(settings)
    teff, logg = output.split('_')[-2:]
    if float(logg) not in settings['C12C13']:
        raise ValueError('C12/C13 ratio for logg={{}} not defined'.format(logg))
    C12C13 = settings['C12C13'][float(logg)]
    rescalc.settings['C12C13'] = C12C13
    rescalc.settings['linelist'] = os.path.realpath(settings['linelist'] + '/C12C13_{{}}'.format(C12C13))
    if not os.path.isdir(rescalc.settings['linelist']):
        raise ValueError('Linelist {{}} not found'.format(rescalc.settings['linelist']))
    if float(logg) not in settings['vturb']:
        raise ValueError('VTURB for logg={{}} not defined'.format(logg))
    rescalc.settings['vturb'] = settings['vturb'][float(logg)]

    output = '{{}}/{{}}'.format(output_dir, output)
    shutil.copytree(output, rundir := (os.path.expandvars('/scratch/$USER/job_$SLURM_JOB_ID/{{}}'.format(os.path.basename(output)))))
    for filename in glob.glob(rescalc.settings['linelist'] + '/*/fort.*') + glob.glob(rescalc.settings['linelist'] + '/*/packager.com') + glob.glob(rescalc.settings['linelist'] + '/*/c12c13.dat'):
        dest = rundir + '/linelist/' + os.path.relpath(filename, rescalc.settings['linelist'])
        os.makedirs(os.path.dirname(dest), exist_ok = True)
        shutil.copy(filename, dest)
    try:
        rescalc.compute_response(rundir + '/response', rundir, rundir + '/linelist')
        shutil.copy('{{}}/response/{{}}_response.pkl'.format(rundir, os.path.basename(output)), '{{}}/response.pkl'.format(output))
        shutil.rmtree(rundir)
        f = open('{{}}/{{}}_complete'.format(results_dir, os.path.basename(output)), 'w'); f.close()
        print('{{}} completed'.format(output), flush = True)
    except Exception as e:
        # shutil.copytree(rundir, output + '/incomplete', symlinks = True)
        shutil.rmtree(rundir)
        f = open('{{}}/{{}}_incomplete'.format(results_dir, os.path.basename(output)), 'w'); f.close()
        raise e
""".strip()


slurm_template = """
#!/bin/bash
#SBATCH --partition={queue}
#SBATCH --account={account}
#SBATCH --job-name="{job_name}"
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={cpus}
#SBATCH --mem={mem}M
#SBATCH -t {walltime}:00:00
#SBATCH --output="{output}"
#SBATCH --error="{error}"
#SBATCH --export=ALL
""".strip()


def dfsynthe():
    tally = {'created': 0, 'successful': 0, 'failed': 0}
    messages = []

    if not os.path.isdir(settings['ODF']['output_dir']):
        os.mkdir(settings['ODF']['output_dir'])

    for population in tqdm.tqdm(griddef.populations):
        # Verify population
        run_dir = '{}/{}/'.format(settings['ODF']['output_dir'], population)
        if os.path.isdir(run_dir):
            validated = os.path.isfile(run_dir + '/kappa.ros')
            try:
                validated = validated and basicatlas.validate_run(run_dir, silent = True)
            except:
                validated = False
            if validated:
                tally['successful'] += 1
            else:
                tally['failed'] += 1
                messages += ['\033[91m{} failed\033[0m'.format(population)]
            continue

        # Calculate new ODF
        f = open('{}/{}.py'.format(settings['scripts_dir'], population), 'w')
        f.write(dfsynthe_template.format(output = run_dir, **griddef.populations[population]))
        f.close()
        tally['created'] += 1

    for message in messages:
        print(message)
    print('\033[91m{failed} failed\033[0m, \033[92m{successful} successful\033[0m, \033[93m{created} created\033[0m of {total}'.format(**tally, total = len(griddef.populations)))

def input_default(prompt, default):
    result = input('{} ({})\n'.format(prompt, default))
    if result == '':
        result = default
    return result

def run():
    # Get job parameters
    params = {
        'size': int(input_default('Batch size?', 1)),
        'ppn': int(input_default('Cores?', 1)),
        'walltime': int(input_default('Walltime?', 24)),
    }

    # Compile the jobs
    script_count = 0
    run_count = 0
    runs = glob.glob('{}/*.py'.format(settings['scripts_dir']))
    for run in runs:
        if run_count == 0:
            script_count += 1
            f = open('{}/run_{}.tqs'.format(settings['runs_dir'], script_count), 'w')
            f.write(slurm_template.format(job_name = 'run_{}'.format(script_count), cpus = params['ppn'], mem = settings['slurm']['mem_per_cpu'] * params['ppn'], walltime = params['walltime'], output = '{}/run_{}.out'.format(settings['runs_dir'], script_count), error = '{}/run_{}.err'.format(settings['runs_dir'], script_count), **settings['slurm']))
            f.write('\n\n')
        run_count += 1
        f.write('python {}\n'.format(run))
        if run_count == params['size']:
            run_count = 0
    if len(runs) > 0:
        f.close()

    # Build the dispatcher
    f = open('{}/launch.com'.format(settings['runs_dir']), 'w')
    runs = glob.glob('{}/*.tqs'.format(settings['runs_dir']))
    for run in runs:
        f.write('sbatch {}\n'.format(run))
    f.close()

def clean():
    files = glob.glob('{}/*'.format(settings['scripts_dir'])) + glob.glob('{}/*'.format(settings['runs_dir']))
    for file in tqdm.tqdm(files):
        os.remove(file)

def clean_results():
    files = glob.glob('{}/*'.format(settings['results_dir']))
    for file in tqdm.tqdm(files):
        os.remove(file)

def reindex_atlas():
    tally = {'created': 0, 'complete': 0, 'deficient': 0, 'excessive': 0}

    available_ATLAS = os.listdir(settings['ATLAS']['output_dir'])
    available_results = os.listdir(settings['results_dir'])

    for population in tqdm.tqdm(griddef.populations):
        expected = [model for model in available_ATLAS if model.find(population) == 0]
        actual = [model for model in available_results if model.find(population) == 0]
        expected_count = len(expected)
        actual_count = len(actual)
        if expected_count == actual_count:
            tally['complete'] += 1
            continue

        if expected_count < actual_count:
            tally['excessive'] += 1
            continue

        tally['created'] += 1
        tally['deficient'] += 1
        f = open('{}/{}.py'.format(settings['scripts_dir'], population), 'w')
        f.write(reindex_atlas_template.format(ATLAS_dir = settings['ATLAS']['output_dir'], results_dir = settings['results_dir'], prefix = population, ATLAS = expected, results = actual))
        f.close()

    print('\033[92m{complete} complete\033[0m, \033[91m{deficient} deficient\033[0m, \033[96m{excessive} excessive\033[0m, \033[93m{created} created\033[0m'.format(**tally))

def atlas():
    tally = {'created': 0, 'converged': 0, 'unconverged': 0, 'failed': 0, 'found': 0}
    messages = []

    if not os.path.isdir(settings['ATLAS']['output_dir']):
        os.mkdir(settings['ATLAS']['output_dir'])

    available_ODFs = os.listdir(settings['ODF']['output_dir'])
    available_ATLAS = set(os.listdir(settings['ATLAS']['output_dir']))
    available_results = os.listdir(settings['results_dir'])
    convergence = {}
    for result in available_results:
        if result[-7:] == '_failed':
            convergence['_'.join(result.split('_')[:-1])] = False
        else:
            convergence['_'.join(result.split('_')[:-2])] = np.array(result.split('_')[-2:]).astype(float)

    for population in tqdm.tqdm(griddef.populations):

        for logg in griddef.logg:


            ##########
            # This is a temporary measure to allow variable carbon sampling
            if logg not in griddef.populations[population]['restrict_logg']:
                continue
            ##########

            # Nominal run parameters
            ODF = '{}/{}/'.format(settings['ODF']['output_dir'], population)
            if population not in available_ODFs:
                messages += ['\033[91m{} misses ODF\033[0m'.format(population)]
                continue
            output = {}; restarts = {}
            teff = np.sort(griddef.teff)[::-1]
            for current, prev in zip(teff, ['auto'] + list(teff[:-1])):
                output[current] = '{}/{}_{}_{}/'.format(settings['ATLAS']['output_dir'], population, current, logg)
                if prev == 'auto':
                    restarts[current] = 'auto'
                else:
                    restarts[current] = '{}/{}_{}_{}/'.format(settings['ATLAS']['output_dir'], population, prev, logg)


            # Evaluate existing structures
            excluded = np.full(len(teff), False)
            for i, model in enumerate(teff):
                filename = output[model]
                while filename[-1] == '/':
                    filename = filename[:-1]
                basename = filename.split('/')[-1]
                if basename in available_ATLAS:
                    excluded[i] = True
                    tally['found'] += 1
                if basename not in convergence and basename in available_ATLAS:
                    messages += ['\033[91m{} did not complete\033[0m'.format(basename)]
                    tally['failed'] += 1
                    continue
                if basename in convergence:
                    if type(convergence[basename]) is bool:
                        messages += ['\033[91mCannot load {}\033[0m'.format(basename)]
                        tally['failed'] += 1
                        continue
                    e1, e2 = convergence[basename]
                    if e1 > settings['ATLAS']['max_flux_error'] or e2 > settings['ATLAS']['max_flux_error_derivative']:
                        # messages += ['\033[91m{} did not converge. Errors {} / {}\033[0m'.format(basename, e1, e2)]
                        tally['unconverged'] += 1
                        continue
                    tally['converged'] += 1
            teff = teff[~excluded]


            # Calculate new structures
            if len(teff) > 0:
                f = open('{}/{}_{}.py'.format(settings['scripts_dir'], population, logg), 'w')
                f.write(atlas_template.format(**griddef.populations[population], restarts = restarts, output = output, logg = logg, teff = '[' + ', '.join(teff.astype(str)) + ']', ODF = ODF, results_dir = settings['results_dir']))
                f.close()
                tally['created'] += len(teff)


    for message in messages:
        print(message)
    print('\033[91m{failed} failed\033[0m, \033[96m{unconverged} unconverged\033[0m, \033[92m{converged} converged\033[0m, \033[93m{created} created\033[0m | \033[94m{found} found\033[0m of {total}'.format(**tally, total = len(griddef.populations) * len(griddef.logg) * len(griddef.teff)))


def delete_atlas():
    delete_hot = True; delete_cold = True; delete_unconv = True
    action = 0
    while True:
        delete_menu = ['Hot failures (Teff > 5000)', 'Cold failures (Teff <= 5000)', 'Include extremely poor convergence', 'Find models', 'Delete models', 'Abort']
        delete_menu[0] = ['[disabled]', '[enabled]'][int(delete_hot)] + ' ' + delete_menu[0]
        delete_menu[1] = ['[disabled]', '[enabled]'][int(delete_cold)] + ' ' + delete_menu[1]
        delete_menu[2] = ['[disabled]', '[enabled]'][int(delete_unconv)] + ' ' + delete_menu[2]
        action = int(pick.pick(delete_menu, indicator = '=>', default_index = action)[1])
        if action == 0:
            delete_hot = not delete_hot
        if action == 1:
            delete_cold = not delete_cold
        if action == 2:
            delete_unconv = not delete_unconv
        if action == 3:
            find_only = True
            break
        if action == 4:
            find_only = False
            break
        if action == 5:
            return

    tally = 0

    available_ATLAS = os.listdir(settings['ATLAS']['output_dir'])
    available_results = os.listdir(settings['results_dir'])
    convergence = {}
    for result in available_results:
        if result[-7:] == '_failed':
            convergence['_'.join(result.split('_')[:-1])] = [result, False]
        else:
            convergence['_'.join(result.split('_')[:-2])] = [result, np.array(result.split('_')[-2:]).astype(float)]

    for model in available_ATLAS:
        teff = float(model.split('_')[-2])
        failed = (model not in convergence) or (type(convergence[model][1]) is bool)
        unconv = (model in convergence) and (type(convergence[model][1]) is not bool) and ((convergence[model][1][0] > 1000) or np.isnan(convergence[model][1][0]))
        if delete_unconv and unconv:
            failed = True
        to_delete = (failed and (teff > 5000) and delete_hot) or (failed and (teff <= 5000) and delete_cold)
        if to_delete:
            tally += 1
            if not find_only:
                shutil.rmtree('{}/{}'.format(settings['ATLAS']['output_dir'], model))
                if model in convergence:
                    os.remove('{}/{}'.format(settings['results_dir'], convergence[model][0]))

    if find_only:
        print('\033[91m{} found\033[0m'.format(tally))
    else:
        print('\033[91m{} deleted\033[0m'.format(tally))

def linelist():
    import rescalc as rescalc_module

    if os.path.isdir(settings['rescalc']['linelist']):
        print('\033[91mExisting linelist found in {}\033[0m'.format(settings['rescalc']['linelist']))
        return
    os.mkdir(settings['rescalc']['linelist'])

    C12C13 = np.unique(list(settings['SYNTHE']['C12C13'].values()))

    for value in C12C13:
        if os.path.isdir(linedir := (settings['rescalc']['linelist'] + '/C12C13_{}'.format(value))):
            continue
        rescalc_module.settings = {**settings['SYNTHE'], **settings['rescalc']}
        rescalc_module.settings['C12C13'] = value
        rescalc_module.settings['vturb'] = 2.0 # Dummy placeholder, VTURB not used in linelist compilation
        rescalc_module.build_all_linelists(linedir)
        print('C12C13={} completed'.format(value))

    return


def rescalc():
    tally = {'successful': 0, 'incomplete': 0, 'excluded': 0, 'created': 0}
    messages = []

    models = set(os.listdir(settings['ATLAS']['output_dir']))
    results = os.listdir(settings['results_dir'])
    results = {'_'.join(result.split('_')[:-1]): result.split('_')[-1] == 'complete' for result in results}

    spectra_per_script = int(input_default('Spectra per script?', 1))
    new = []

    # We are computing models out of order to ensure that the sampling is as uniform as possible when incomplete
    shuffled_models = list(models)
    np.random.shuffle(shuffled_models)
    for model in tqdm.tqdm(shuffled_models):
        if model in results and results[model]:
            tally['successful'] += 1
            continue

        elif model in results and (not results[model]):
            # messages += ['\033[91m{} is incomplete\033[0m'.format(model)]
            tally['incomplete'] += 1
            continue

        teff, logg = model.split('_')[-2:]

        # This is for excluding certain models from TiO calculation
        # zscale = float(model.split('_')[0][1:])
        # alpha = float(model.split('_')[1][1:])
        # if int(teff) > 6000 or (int(teff) > 5000 and zscale + alpha < 0):
        #     tally['excluded'] += 1
        #     continue

        new += [model]

    # Calculate new spectra
    new = [new[i:i + spectra_per_script] for i in range(0, len(new), spectra_per_script)]
    rescalc_settings = {**settings['SYNTHE'], **settings['rescalc']}
    rescalc_path = os.path.dirname(os.path.realpath(__file__))
    for i, script in enumerate(tqdm.tqdm(new)):
        f = open('{}/script_{}.py'.format(settings['scripts_dir'], i), 'w')
        f.write(rescalc_template.format(output = script, settings = rescalc_settings, rescalc_path = rescalc_path, output_dir = settings['ATLAS']['output_dir'], results_dir = settings['results_dir']))
        f.close()
        tally['created'] += len(script)

    for message in messages:
        print(message)
    print('\033[92m{successful} successful\033[0m, \033[96m{excluded} excluded\033[0m, \033[91m{incomplete} incomplete\033[0m, \033[93m{created} created\033[0m of {total}'.format(**tally, total = len(models)))







#############################################################################################


globals()[menu[action][1]]()
