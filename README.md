# orderly

This is an [`orderly`](https://github.com/vimc/orderly) project.  The directories are:

* `src`: create new reports here
* `archive`: versioed results of running your report
* `data`: copies of data used in the reports

(you can delete or edit this file safely)

## Instructions for DIDE server

#### 1. pars_init generation




## Instructions for Azure server [OLD]

This is still a manual process requiring logging into the server and executing the model fits.

Running everything takes a while so you really probably want to run this under screen. After logging in, you can start a screen session with:

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

## Docker Images and Pulling New Code Changes

If there are code changes, then the images will take a while to build. These are separated into two components "base", which takes quite a while to build, and another image built on top of this that includes the fast moving packages.

These are built automatically on [buildkite](https://buildkite.com/mrc-ide/global-lmic-report) but will take a number of minutes to complete.  You can do this manually with:

```
./docker/build-base
./docker/build
```

Typically the second will be enough unless new packages have been added to the base dockerfile, and that step runs in a few seconds.

## Running reports

Run all the reports and build index pages with

```
./docker/run_google_pmcmc_spline_np

```

If you want/need to run a previous date then pass that as an argument

```
./docker/run_google_pmcmc_spline_np 2020-11-10
```

If a particular country run fails (can happen due to data streams changing) and 
you need to rerun a specific country :

```
./docker/run_country_google_pmcmc_spline_np 2020-11-10 UZB
```

If you need to run all the countries again for a particular date (e.g. squire has
updated and the drat repo (https://ncov-ic.github.io/drat) was not updated):


```
./docker/run_google_countries_pmcmc_spline_np 2020-11-10
```

This will run all the model fits, but won't collate the website together. This 
can be achieved using:

```
./docker/run_collate_google_pmcmc_spline_np 2020-11-10
```

If you need to only run some of the countries again for a particular date:

```
echo -en "UZB\nGBR" > redo
./docker/run_google_countries_pmcmc_spline_np 2020-11-10 FALSE FALSE TRUE TRUE redo
```

There are 4 other arguments in the run functions (as shown above) The default 
arguments are shown above. See specific `run_<>.sh` files for more details.

## Publishing reports

To publish the website, you will need the ssh key (this has been done on the server already)

```
./scripts/get_deploy/key
```

Deploy to staging with (n.b. need to add in the staging navbar code)

```
./scripts/publish_website staging
```

Check the [staging website](https://mrc-ide.github.io/global-lmic-reports-staging/)

Deploy to production with

```
./scripts/publish_website production
```
