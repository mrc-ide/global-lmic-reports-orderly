# orderly

This is an [`orderly`](https://github.com/vimc/orderly) project.  The directories are:

* `src`: create new reports here
* `archive`: versioed results of running your report
* `data`: copies of data used in the reports

(you can delete or edit this file safely)

## Instructions

Currently somewhat manual, we'll move to an automated system once the data sources fully settle down.

If there are code changes, then the images will take a while to build. These are separated into two components "base", which takes quite a while to build, and another image built on top of this that includes the fast moving packages.

```
./docker/build-base
./docker/build
```

Run all the reports and build index pages

```
./docker/run
```

If you want/need to run a previous date then pass that as an argument

```
./docker/run 2020-04-01
```

To publish the website, you will need the ssh key (this has been done on the server already)

```
./scripts/get_deploy/key
```

Deploy to staging with

```
./scripts/publish_website staging
```

Deploy to production with

```
./scripts/publish_website production
```
