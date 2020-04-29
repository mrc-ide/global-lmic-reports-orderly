# orderly

This is an [`orderly`](https://github.com/vimc/orderly) project.  The directories are:

* `src`: create new reports here
* `archive`: versioed results of running your report
* `data`: copies of data used in the reports

(you can delete or edit this file safely)

## Instructions on the server

Currently somewhat manual, we'll move to an automated system once the data sources fully settle down.

Running everything takes a while so you really probably want to run this under screen:

```
screen
```

if your network connection is disrupted then things will continue.  You can connect to an existing screen using:

```
screen -r
```

You will probably end up with a surplus of old screen commands around, which you can see with

```
screen -list
```

These commands might be useful to clean them out:

```
screen -X -S [session # you want to kill] quit
screen -wipe
```

If there are code changes, then the images will take a while to build. These are separated into two components "base", which takes quite a while to build, and another image built on top of this that includes the fast moving packages.

These are built automatically on [buildkite](https://buildkite.com/mrc-ide/global-lmic-report) but will take a number of minutes to complete.  You can do this manually with:

```
./docker/build-base
./docker/build
```

Typically the second will be enough unless new packages have been added to the base dockerfile, and that step runs in a few seconds.

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

Check the [staging website](https://mrc-ide.github.io/global-lmic-reports-staging/)

Deploy to production with

```
./scripts/publish_website production
```
