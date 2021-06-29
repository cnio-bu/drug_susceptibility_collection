import pandas as pd

datasets = pd.read_csv('datasets.csv', sep=',', index_col=0)
results  = '/raid/sagarcia/test_drugs'


def get_prism_genesets(wildcards):
    checkpoint_output = checkpoints.prism_annotate_models.get(**wildcards).output['auc_models_candidates']

    return expand(f'{results}/prism/genesets/{{broad_id}}',
           broad_id=glob_wildcards(os.path.join(checkpoint_output, "{brd_id}.csv")).brd_id)


def get_gdsc_ebayes(wildcards):
    checkpoint_output = checkpoints.gdsc_generate_compound_curves.get(**wildcards).output['auc_models_candidates']

    return expand(f'{results}/gdsc/ebayes/{{drug_id}}_eBayes.rds',
           drug_id=glob_wildcards(os.path.join(checkpoint_output, "{drug_id}.csv")).drug_id)


rule all:
    input:
       # get_prism_genesets,
       get_gdsc_ebayes




## Load rules ##
include: 'rules/common.smk'
include: 'rules/gdsc_arrays_signatures.smk'
include: 'rules/prism_signatures.smk'

