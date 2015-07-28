web: gunicorn app:app --log-file=-
celeryd: celery -A tasks worker --loglevel=info -E