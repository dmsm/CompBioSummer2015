# CompBioSummer2015
Research in Computational Biology at Harvey Mudd College, Summer 2015

## Installing on Heroku
First create a Heroku app with `heroku create` in the directory of the repo.

Set up an Amazon S3 bucket, and configure Heroku config vars `AWS_ACCESS_KEY`, `AWS_SECRET_ACCESS_KEY`, and `AWS_SECRET_ACCESS_KEY` accordingly with
```
heroku config:set VAR_NAME=var_value
```

Specify a custom buildpack by running
```
heroku buildpacks:set git://github.com/dulaccc/heroku-buildpack-geodjango.git#1.1
```
and provision a Redis instance with a
```
heroku addons:create redistogo
```

Now, run `git push heroku master` and scale the app and worker:
```
heroku ps:scale web=1
heroku scale worker=1
```