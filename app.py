#!/usr/local/bin/python
import os

from flask import Flask, request, render_template, send_from_directory
from werkzeug.utils import secure_filename
from rq import Queue
import tinys3

from tasks import process_files
from worker import conn

app = Flask(__name__)
q = Queue(connection=conn)

UPLOAD_FOLDER = "/app/temp"
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
        # handle uploaded file
        newick_file = request.files['newick']
        if newick_file and allowed_file(newick_file.filename):
            conn = tinys3.Connection(os.environ.get('AWS_ACCESS_KEY'), os.environ.get('AWS_SECRET_ACCESS_KEY'),
                                     default_bucket=os.environ.get('S3_BUCKET_NAME'))
            filename = secure_filename(newick_file.filename)
            conn.upload(filename, newick_file)
        else:
            return render_template("documentation.html")

        dup = request.form['dup'] if request.form['dup'] != '' else 2
        trans = request.form['trans'] if request.form['dup'] != '' else 3
        loss = request.form['loss'] if request.form['dup'] != '' else 1

        switch_hi = request.form['switchhigh'] if request.form['dup'] != '' else 4.5
        switch_lo = request.form['switchlow'] if request.form['dup'] != '' else 1.5

        loss_hi = request.form['losshigh'] if request.form['dup'] != '' else 3
        loss_lo = request.form['losslow'] if request.form['dup'] != '' else 1

        job = q.enqueue_call(func=process_files,
                             args=(filename, dup, trans, loss, request.form['scoring'], switch_lo, switch_hi, loss_lo,
                                   loss_hi),
                             timeout=600)

        return render_template("results.html", task_id=job.id)


@app.route('/status/<task_id>')
def taskstatus(task_id):
    job = q.fetch_job(task_id)
    if job.result:
        return render_template("display.html", **job.result)
    elif:
        job.
    else:
        return "PENDING"


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