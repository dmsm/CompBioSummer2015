#!/usr/local/bin/python
import os
import ast

from flask import Flask, request, render_template, send_from_directory
from werkzeug.utils import secure_filename
from compbio.vis import transsvg
from rasmus import treelib1
from compbio import phylo
from celery import Celery

from MasterReconciliation import Reconcile
from ReconConversion import freqSummation

app = Flask(__name__)

app.config['BROKER_URL'] = os.environ['REDIS_URL']
app.config['CELERY_RESULT_BACKEND'] = os.environ['REDIS_URL']
celery = Celery(app.name, broker=app.config['REDIS_URL'])
celery.conf.update(app.config)

UPLOAD_FOLDER = "tmp"
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER


@app.route('/')
def index():
    """ Returns the index page"""
    return render_template('index.html', page='home')


@app.route('/form')
def form():
    """Returns the form page"""
    return render_template('form.html', page='upload')


@app.route('/documentation')
def documentation():
    """Returns the documentation page"""
    return render_template('documentation.html')


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] == "newick"


@celery.task(bind=True)
def process_files(self, dup, trans, loss, *args):
    raw_name = os.path.splitext(os.path.basename(args[0]))[0]
    Reconcile(args)
    freqSummation(args)
    with open(os.path.join(app.config['UPLOAD_FOLDER'], "{}freqFile.txt".format(raw_name))) as f:
        lines = f.readlines()

    score_list = ast.literal_eval(lines[0])
    total_freq = float(lines[1][:-2])
    total_recon = float(lines[3])
    total_cost = float(lines[2][:-2])

    if request.form['scoring'] == "Frequency":
        score_method = "Frequency"
    elif request.form['scoring'] == "xscape":
        score_method = "Xscape Scoring"
    else:
        score_method = "Unit Scoring"

    results_list = []
    for x, score in enumerate(score_list):
        tree = treelib1.read_tree(os.path.join(app.config['UPLOAD_FOLDER'], "{}.tree".format(raw_name)))
        stree = treelib1.read_tree(os.path.join(app.config['UPLOAD_FOLDER'], "{}{}.stree".format(raw_name, x)))
        brecon = phylo.read_brecon(os.path.join(app.config['UPLOAD_FOLDER'],
                                                "{}{}.mowgli.brecon".format(raw_name, x)), tree, stree)
        output = os.path.join(app.config['UPLOAD_FOLDER'], "{}{}.svg".format(raw_name, x))
        phylo.add_implied_spec_nodes_brecon(tree, brecon)
        transsvg.draw_tree(tree, brecon, stree, filename=output)

        percent = 100.0 * score / total_freq
        running_tot_score = sum(score_list[:x])
        running_tot = min(100.0 * running_tot_score / total_freq, 100)

        results_list.append((score, percent, running_tot))

    self.update_state(state='SUCCESS')

    return {results_list: results_list, raw_name: raw_name, dup: dup, trans: trans, loss: loss, total_cost: total_cost,
            score_method: score_method, total_freq: total_freq, total_recon: total_recon}



@app.route('/reconcile', methods=['POST'])
def reconcile(carousel=None):
    """ Creates the results page using MasterReconciliation and vistrans"""
    if request.method == 'POST':
        # create upload dir if doesn't exits
        if not os.path.exists(app.config['UPLOAD_FOLDER']):
            os.makedirs(app.config['UPLOAD_FOLDER'])

        # clear out files from last run
        files = os.listdir(app.config['UPLOAD_FOLDER'])
        for f in files:
            os.remove(os.path.join(app.config['UPLOAD_FOLDER'], f))

        # handle uploaded file
        newick_file = request.files['newick']
        if newick_file and allowed_file(newick_file.filename):
            filename = secure_filename(newick_file.filename)
            # raw_name = os.path.splitext(os.path.basename(filename))[0]
            file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            newick_file.save(file_path)
        else:
            return render_template("documentation.html")

        dup = request.form['dup'] if request.form['dup'] != '' else 2
        trans = request.form['trans'] if request.form['dup'] != '' else 3
        loss = request.form['loss'] if request.form['dup'] != '' else 1

        switch_hi = request.form['switchhigh'] if request.form['dup'] != '' else 4.5
        switch_lo = request.form['switchlow'] if request.form['dup'] != '' else 1.5

        loss_hi = request.form['losshigh'] if request.form['dup'] != '' else 3
        loss_lo = request.form['losslow'] if request.form['dup'] != '' else 1

        task = process_files.apply_async(args=[dup, trans, loss, file_path, dup, trans, loss, request.form['scoring'], switch_lo, switch_hi, loss_lo, loss_hi])
        # Reconcile([file_path, dup, trans, loss, request.form['scoring'], switch_lo, switch_hi, loss_lo, loss_hi])
        # freqSummation([file_path, dup, trans, loss, request.form['scoring'], switch_lo, switch_hi, loss_lo, loss_hi])

        # with open(os.path.join(app.config['UPLOAD_FOLDER'], "{}freqFile.txt".format(raw_name))) as f:
        #     lines = f.readlines()
        #
        # score_list = ast.literal_eval(lines[0])
        # total_freq = float(lines[1][:-2])
        # total_recon = float(lines[3])
        # total_cost = float(lines[2][:-2])
        #
        # if request.form['scoring'] == "Frequency":
        #     score_method = "Frequency"
        # elif request.form['scoring'] == "xscape":
        #     score_method = "Xscape Scoring"
        # else:
        #     score_method = "Unit Scoring"
        #
        # results_list = []
        # for x, score in enumerate(score_list):
        #     tree = treelib1.read_tree(os.path.join(app.config['UPLOAD_FOLDER'], "{}.tree".format(raw_name)))
        #     stree = treelib1.read_tree(os.path.join(app.config['UPLOAD_FOLDER'], "{}{}.stree".format(raw_name, x)))
        #     brecon = phylo.read_brecon(os.path.join(app.config['UPLOAD_FOLDER'],
        #                                             "{}{}.mowgli.brecon".format(raw_name, x)), tree, stree)
        #     output = os.path.join(app.config['UPLOAD_FOLDER'], "{}{}.svg".format(raw_name, x))
        #     phylo.add_implied_spec_nodes_brecon(tree, brecon)
        #     transsvg.draw_tree(tree, brecon, stree, filename=output)
        #
        #     percent = 100.0 * score / total_freq
        #     running_tot_score = sum(score_list[:x])
        #     running_tot = min(100.0 * running_tot_score / total_freq, 100)
        #
        #     results_list.append((score, percent, running_tot))

        # return render_template("results.html", results_list=results_list, raw_name=raw_name,
        #                        dup=dup, trans=trans, loss=loss, total_cost=total_cost, score_method=score_method,
        #                        total_freq=total_freq, total_recon=total_recon)

        return render_template("results.html", task_id=task.id)



@app.route('/status/<task_id>')
def taskstatus(task_id):
    task = process_files.AsyncResult(task_id)
    if task.state == 'PENDING':
        return "PENDING"
    else:
        return render_template("display.html", **task.info)


@app.route('/uploads/<filename>')
def send_file(filename):
    """Takes in a filename and sends it from the directory to the results page"""
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename)


@app.after_request
def add_header(response):
    """
    Add headers to both force latest IE rendering engine or Chrome Frame,
    and also to cache the rendered page for 10 minutes.
    """
    response.headers['X-UA-Compatible'] = 'IE=Edge,chrome=1'
    response.headers['Cache-Control'] = 'public, max-age=0'
    return response

app.debug = True
if __name__ == '__main__':
    app.run()