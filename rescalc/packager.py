from config import settings
import griddef
import numpy as np
import tqdm
import os

import pick

# Main menu options
menu = [
    ['Package ODFs', 'dfsynthe'],
    ['Package structures and spectra', 'atlas'],
]

# What are we doing?
action = int(pick.pick([menu[i][0] for i in range(len(menu))], indicator = '=>')[1])

template_dfsynthe = """
import os, shutil

input = '{input}'
output = '{output}'
staging = os.path.expandvars('/scratch/$USER/job_$SLURM_JOB_ID/{{}}'.format(os.path.basename(input)))

if os.path.isdir(staging):
    shutil.rmtree(staging)
os.mkdir(staging)

for file in ['p00big2.bdf', 'kappa.ros', 'xnfdf.com', 'xnfdf.out']:
    shutil.copy('{{}}/{{}}'.format(input, file), staging)
    assert os.path.getsize('{{}}/{{}}'.format(staging, file)) > 0

f = open(filename := ('{{}}/{{}}'.format(staging, 'xnfdf.out')), 'r')
content = f.read()
f.close()
f = open(filename, 'w')
start = content.rfind('\\n', 0, content.rfind('\\n', 0, content.rfind('\\n', 0, content.find('0TITLE'))))
end = content.find('\\n', content.find('\\n', content.find('\\n', content.find('0FREQID'))))
f.write(content[start:end])
f.close()

os.system('cd {{staging}}/.. && zip -r {{output}} {{ODF}}/*'.format(output = output, staging = staging, ODF = os.path.basename(staging)))
print(os.path.basename(input) + ' complete')

""".strip()

template_atlas = """

import os, shutil, pickle
import numpy as np

model_dir = '{model_dir}'
population = '{population}'
models = {models}

workdir = os.path.expandvars('/scratch/$USER/job_$SLURM_JOB_ID/{{}}'.format(population))
if os.path.isdir(workdir):
    shutil.rmtree(workdir)
os.makedirs(structures_dir := (workdir + '/structures/'), exist_ok = True)
os.makedirs(fullres_dir := (workdir + '/fullres/{{}}/'.format(population)), exist_ok = True)
os.makedirs(binned_dir := (workdir + '/binned/{{}}/'.format(population)), exist_ok = True)

total_size = 0
binned_size = 0

for model in models:
    origin = '{{}}/{{}}/response.pkl'.format(model_dir, model)
    if not os.path.isfile(origin):
        continue
    shutil.copy(origin, destination := ('{{}}/{{}}.pkl'.format(fullres_dir, model)))
    f = open(destination, 'rb')
    response = pickle.load(f)
    f.close()
    wl = np.exp(np.arange(np.ceil(np.log(response['meta']['wl_start']) / (lgr := np.log(1.0 + 1.0 / response['meta']['res']))), np.floor(np.log(response['meta']['wl_end']) / lgr) + 1) * lgr) * 10
    assert len(response['null']['line']) == len(wl)
    assert len(response['null']['cont']) == len(wl)
    assert not np.any(np.isnan(response['null']['line']))
    assert not np.any(np.isnan(response['null']['cont']))
    assert np.min(response['null']['cont']) >= 0
    assert np.min(response['null']['line']) >= 0
    total_size += os.path.getsize(destination)

    # Produce a binned version of the spectrum
    binning = 10
    response['meta']['binning'] = binning
    response['null']['line'] = response['null']['line'][:(len(response['null']['line']) // binning) * binning].reshape(-1, binning).mean(axis = 1)
    response['null']['cont'] = response['null']['cont'][:(len(response['null']['cont']) // binning) * binning].reshape(-1, binning).mean(axis = 1)
    f = open(destination := ('{{}}/{{}}.pkl'.format(binned_dir, model)), 'wb')
    pickle.dump(response, f)
    f.close()
    binned_size += os.path.getsize(destination)

    # Save the structure in the BasicATLAS format
    os.makedirs(destination := ('{{}}/{{}}/'.format(structures_dir, model)))
    shutil.copy('{{}}/{{}}/output_summary.out'.format(model_dir, model), destination)
    shutil.copy('{{}}/{{}}/output_last_iteration.out'.format(model_dir, model), destination)
    f = open('{{}}/output_main.out'.format(destination), 'w'); f.close()

output_dir = os.path.realpath('{{}}/../output/fullres/'.format(model_dir))
os.makedirs(output_dir, exist_ok = True)
f = open(filename := (fullres_dir + '/../compressor.com'), 'w')
f.write('cd {{}}\\n'.format(os.path.realpath(fullres_dir + '/..')))
f.write('zip -q -r {{output_dir}}/{{population}}.zip {{population}}\\n'.format(output_dir = output_dir, population = population))
f.write('zip -T {{output_dir}}/{{population}}.zip > compressor.out\\n'.format(output_dir = output_dir, population = population))
f.write('md5sum {{output_dir}}/{{population}}.zip >> compressor.out\\n'.format(output_dir = output_dir, population = population))
f.close()
os.system('bash {{}}'.format(filename))
f = open(fullres_dir + '/../compressor.out', 'r')
result = f.read().strip()
f.close()
if result.find(' OK') == -1:
    raise ValueError('ZIP corrupted')
print(result)
print('Fullres filesize: {{}}'.format(total_size))

output_dir = os.path.realpath('{{}}/../output/binned/'.format(model_dir))
os.makedirs(output_dir, exist_ok = True)
f = open(filename := (binned_dir + '/../compressor.com'), 'w')
f.write('cd {{}}\\n'.format(os.path.realpath(binned_dir + '/..')))
f.write('zip -q -r {{output_dir}}/{{population}}.zip {{population}}\\n'.format(output_dir = output_dir, population = population))
f.write('zip -T {{output_dir}}/{{population}}.zip > compressor.out\\n'.format(output_dir = output_dir, population = population))
f.write('md5sum {{output_dir}}/{{population}}.zip >> compressor.out\\n'.format(output_dir = output_dir, population = population))
f.close()
os.system('bash {{}}'.format(filename))
f = open(binned_dir + '/../compressor.out', 'r')
result = f.read().strip()
f.close()
if result.find(' OK') == -1:
    raise ValueError('ZIP corrupted')
print(result)
print('Binned filesize: {{}}'.format(binned_size))

output_dir = os.path.realpath('{{}}/../output/structures/'.format(model_dir))
os.makedirs(output_dir, exist_ok = True)
f = open(filename := (structures_dir + '/compressor.com'), 'w')
f.write('cd {{}}\\n'.format(os.path.realpath(structures_dir)))
f.write('zip -q -r {{output_dir}}/{{population}}.zip z*\\n'.format(output_dir = output_dir, population = population))
f.write('zip -T {{output_dir}}/{{population}}.zip > compressor.out\\n'.format(output_dir = output_dir, population = population))
f.write('md5sum {{output_dir}}/{{population}}.zip >> compressor.out\\n'.format(output_dir = output_dir, population = population))
f.close()
os.system('bash {{}}'.format(filename))
f = open(structures_dir + '/compressor.out', 'r')
result = f.read().strip()
f.close()
if result.find(' OK') == -1:
    raise ValueError('ZIP corrupted')
print(result)

shutil.rmtree(workdir)

""".strip()


def dfsynthe():
    ODFs = os.listdir(settings['ODF']['output_dir'])
    for ODF in tqdm.tqdm(ODFs):
        output = '/expanse/lustre/projects/csd835/rgerasim/pfsgrid/temp/{}.zip'.format(ODF)
        f = open('/expanse/lustre/projects/csd835/rgerasim/pfsgrid/temp/script.py', 'w')
        f.write(template_dfsynthe.format(input = settings['ODF']['output_dir'] + '/' + ODF, output = output))
        f.close()
        break

def atlas():
    models = os.listdir(settings['ATLAS']['output_dir'])
    if os.path.isdir(existing := (settings['ATLAS']['output_dir'] + '/../output/fullres')):
        existing = [model[:-4] for model in os.listdir(existing)]
    else:
        existing = []

    populations = {}

    Cgrid = {}

    for model in tqdm.tqdm(models):
        components = model.split('_')
        Cgrid_id = components[0] + '_' + components[-1]
        if Cgrid_id not in Cgrid:
            Cgrid[Cgrid_id] = np.round([griddef.interpolator([float(components[0][1:]), float(components[-1]), i])[0] + 0.0001 for i in np.arange(-2, 3)], 1)
        carbon = np.where(np.unique(Cgrid[Cgrid_id]) == float(components[2][1:]))[0][0]
        population = '{}_{}_c{}'.format(components[0], components[1], ['ud', 'd', 'm', 'e', 'ue'][carbon])
        if population not in populations:
            populations[population] = []
        populations[population] += [model]

    assert len(populations) == len(griddef.zscale) * len(griddef.alpha) * len(list(Cgrid.values())[0])
    for population in populations:
        assert len(populations[population]) == len(np.unique(populations[population]))
        assert len(populations[population]) == len(griddef.logg) * len(griddef.teff)

    for population in tqdm.tqdm(populations):
        if population in existing:
            continue
        f = open('{}/{}.py'.format(settings['scripts_dir'], population), 'w')
        f.write(template_atlas.format(model_dir = settings['ATLAS']['output_dir'], models = populations[population], population = population))
        f.close()

globals()[menu[action][1]]()
