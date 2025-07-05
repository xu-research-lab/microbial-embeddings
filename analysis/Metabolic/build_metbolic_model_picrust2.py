from carveme import config, project_dir
from carveme import __version__ as version
from carveme.reconstruction.carving import carve_model, build_ensemble
from carveme.reconstruction.eggnog import load_eggnog_data
from carveme.reconstruction.gapfilling import multiGapFill
from carveme.reconstruction.utils import load_media_db, load_soft_constraints, load_hard_constraints, annotate_genes
from carveme.reconstruction.ncbi_download import load_ncbi_table, download_ncbi_genome
from carveme.reconstruction.scoring import reaction_scoring
from reframed.cobra.ensemble import save_ensemble
from reframed import load_cbmodel, save_cbmodel, Environment, set_default_solver
from reframed.io.sbml import sanitize_id
from reframed.core.transformation import apply_bounds
import argparse
import os
import os.path
import pandas as pd
from multiprocessing import Pool
from glob import glob
import subprocess
import sys
import json

from joblib import Parallel, delayed
from tqdm import tqdm
import time

def maincall(inputfile, outputfile, model_id, default_score=-1.0, uptake_score=0.0, soft_score=1.0,
             ref_score=0.0, universe=None, universe_file=None, ensemble_size=None, 
             verbose=False, debug=False, flavor=None, gapfill=None, blind_gapfill=False, init=None,
             mediadb=None, soft=None, hard=None, reference=None, recursive_mode=False):

    if soft:
        try:
            soft_constraints = load_soft_constraints(soft)
        except IOError:
            raise IOError('Failed to load soft-constraints file:' + soft)
    else:
        soft_constraints = None

    if hard:
        try:
            hard_constraints = load_hard_constraints(hard)
        except IOError:
            raise IOError('Failed to load hard-constraints file:' + hard)
    else:
        hard_constraints = None

    if reference:
        if verbose:
            print('Loading reference model...')

        try:
            ref_model = load_cbmodel(reference)
        except:
            raise IOError('Failed to load reference model.')
    else:
        ref_model = None
    

    annotations = pd.read_csv(inputfile)

    if not universe_file:
        if universe:
            universe_file = f"{project_dir}{config.get('generated', 'folder')}universe_{universe}.xml.gz"
        else:
            universe_file = project_dir + config.get('generated', 'default_universe')

    try:
        universe_model = load_cbmodel(universe_file, flavor='bigg')
        universe_model.id = model_id
        apply_bounds(universe_model)
    except IOError:
        available = '\n'.join(glob(f"{config.get('generated', 'folder')}universe_*.xml.gz"))
        raise IOError(f'Failed to load universe model: {universe_file}\nAvailable universe files:\n{available}')

    if gapfill or init:

        if verbose:
            print('Loading media library...')

        if not mediadb:
            mediadb = project_dir + config.get('input', 'media_library')

        try:
            media_db = load_media_db(mediadb)
        except IOError:
            raise IOError('Failed to load media library:' + mediadb)

    if verbose:
        print('Scoring reactions...')

    gene_annotations = pd.read_csv(project_dir + config.get('generated', 'gene_annotations'), sep='\t')
    bigg_gprs = project_dir + config.get('generated', 'bigg_gprs')
    gprs = pd.read_csv(bigg_gprs)
    gprs = gprs[gprs.reaction.isin(universe_model.reactions)]

    debug_output = model_id if debug else None
    scores, gene2gene = reaction_scoring(annotations, gprs, debug_output=debug_output)

    if scores is None:
        print('The input genome did not match sufficient genes/reactions in the database.')
        return

    if not flavor:
        flavor = config.get('sbml', 'default_flavor')

    init_env = None

    if init:
        if init in media_db:
            init_env = Environment.from_compounds(media_db[init])
        else:
            print(f'Error: medium {init} not in media database.')

    universe_model.metadata['Description'] = 'This model was built with CarveMe version ' + version

    if ensemble_size is None or ensemble_size <= 1:
        if verbose:
            print('Reconstructing a single model')

        model = carve_model(universe_model, scores, inplace=(not gapfill), default_score=default_score,
                            uptake_score=uptake_score, soft_score=soft_score, soft_constraints=soft_constraints,
                            hard_constraints=hard_constraints, ref_model=ref_model, ref_score=ref_score,
                            init_env=init_env, debug_output=debug_output)
        annotate_genes(model, gene2gene, gene_annotations)

    else:
        if verbose:
            print('Building an ensemble of', ensemble_size, 'models')

        ensemble = build_ensemble(universe_model, scores, ensemble_size, init_env=init_env)

        annotate_genes(ensemble.model, gene2gene, gene_annotations)
        save_ensemble(ensemble, outputfile, flavor=flavor)
        return

    if model is None:
        print("Failed to build model.")
        return

    if not gapfill:
        save_cbmodel(model, outputfile, flavor=flavor)

    else:
        media = gapfill.split(',')

        if verbose:
            m1, n1 = len(model.metabolites), len(model.reactions)
            print(f"Gap filling for {', '.join(media)}...")

        max_uptake = config.getint('gapfill', 'max_uptake')

        if blind_gapfill:
            scores = None
        else:
            scores = dict(scores[['reaction', 'normalized_score']].values)
        multiGapFill(model, universe_model, media, media_db, scores=scores, max_uptake=max_uptake, inplace=True)

        if verbose:
            m2, n2 = len(model.metabolites), len(model.reactions)
            print(f'Added {(n2 - n1)} reactions and {(m2 - m1)} metabolites')

        if init_env:  # Initializes environment again as new exchange reactions can be acquired during gap-filling
            init_env.apply(model, inplace=True, warning=False)

        save_cbmodel(model, outputfile, flavor=flavor)

    if verbose:
        print('Done.')


def main(i):
    model_id = i
    inputfile = f"Data/OTU_bigg_gene/{model_id}.tsv"
    outputfile = f"Data/OTU_metabolic_model_M3/{model_id}.xml"
    universe = fid_gram.get(model_id)
    mediadb = "Data/media_db.tsv"
    maincall(inputfile=inputfile, outputfile=outputfile, model_id=model_id,
                universe=universe, gapfill="M3", mediadb=mediadb)


if __name__ == '__main__':
    set_default_solver("cplex")
    with open("Data/fid_gram.json", "r") as file:
        fid_gram = json.load(file) 
    fid = sys.argv[1]
    main(fid)

