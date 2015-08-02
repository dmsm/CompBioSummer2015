import os
import urllib
import ast

import tinys3
from compbio.vis import transsvg
from rasmus import treelib1
from compbio import phylo

from MasterReconciliation import Reconcile
from ReconConversion import freqSummation


def process_files(*args):
    conn = tinys3.Connection(os.environ.get('AWS_ACCESS_KEY'), os.environ.get('AWS_SECRET_ACCESS_KEY'),
                             default_bucket=os.environ.get('S3_BUCKET_NAME'))
    filename = args[0]
    url = "http://s3.amazonaws.com/{}/{}".format(os.environ.get('S3_BUCKET_NAME'), filename)
    urllib.urlretrieve(url, filename)
    raw_name = os.path.splitext(filename)[0]
    Reconcile(args)
    freqSummation(args)
    with open("{}freqFile.txt".format(raw_name)) as f:
        lines = f.readlines()

    score_list = ast.literal_eval(lines[0])
    total_freq = float(lines[1][:-2])
    total_recon = float(lines[3])
    total_cost = float(lines[2][:-2])

    scoring = args[4]
    if scoring == "Frequency":
        score_method = "Frequency"
    elif scoring == "xscape":
        score_method = "Xscape Scoring"
    else:
        score_method = "Unit Scoring"

    results_list = []
    for x, score in enumerate(score_list):
        tree = treelib1.read_tree("{}.tree".format(raw_name))
        stree = treelib1.read_tree("{}{}.stree".format(raw_name, x))
        brecon = phylo.read_brecon("{}{}.mowgli.brecon".format(raw_name, x), tree, stree)
        output = "/{}{}.svg".format(raw_name, x)
        phylo.add_implied_spec_nodes_brecon(tree, brecon)
        transsvg.draw_tree(tree, brecon, stree, filename=output)
        f = open(output, 'rb')
        conn.upload(filename, f)

        percent = 100.0 * score / total_freq
        running_tot_score = sum(score_list[:x])
        running_tot = min(100.0 * running_tot_score / total_freq, 100)

        results_list.append((score, percent, running_tot))

    return {'results_list': results_list, 'raw_name': raw_name, 'dup': args[1], 'trans': args[2], 'loss': args[3],
            'total_cost': total_cost, 'score_method': score_method, 'total_freq': total_freq,
            'total_recon': total_recon}

