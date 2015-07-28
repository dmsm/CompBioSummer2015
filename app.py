#!/usr/local/bin/python
import os

from flask import Flask, request, render_template, send_from_directory
from werkzeug.utils import secure_filename
from rq import Queue

from tasks import process_files
from worker import conn

app = Flask(__name__)
q = Queue(connection=conn)

UPLOAD_FOLDER = "/tmp"
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

        job = q.enqueue(process_files, app.config['UPLOAD_FOLDER'], dup, trans, loss, request.form['scoring'], file_path, dup, trans, loss, request.form['scoring'], switch_lo, switch_hi, loss_lo, loss_hi)

        return "{} {} {}".format(job.id, file_path, open(file_path, 'r').name)



@app.route('/status/<task_id>')
def taskstatus(task_id):
    job = q.fetch_job(task_id)
    if job.result:
        return render_template("display.html", **job.result)
    else:
        return "PENDING"


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