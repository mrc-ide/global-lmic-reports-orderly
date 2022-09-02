# orderly

This is an [`orderly`](https://github.com/vimc/orderly) project.  The directories are:

* `src`: create new reports here
* `archive`: versioned results of running your report
* `vignettes`: `.Rmd` files for running the fitting process 
* `R`: Functions used in the vignettes
* `global`: Resources used for setting up the GitHub page, accessible to orderly tasks
* `scripts`: `.sh` scripts for updating the output GitHub repositories

The following directories will be generated during regular usage:

* `fits`: `.pdf` files for checking model fit quality.
* `draft`: draft versions of the orderly tasks.
* `gh-esft`: copy of the output repository [global_lmic_projections_esft](https://github.com/mrc-ide/global_lmic_projections_esft)
* `gh-fits`: copy of the output repository [nimue_global_fits](https://github.com/mrc-ide/nimue_global_fits)
* `gh-pages`: copy of the output repository [global-lmic-reports](https://github.com/mrc-ide/global-lmic-reports)

# Instructions for DIDE server

On overview of the submission of model fits to the DIDE cluster is in the 
`dide_run.Rmd` vignette in the vignettes folder. All instructions should be self
contained in that Rmd document, and it is advised to follow through that document. 
To see an overview, click [here](https://htmlpreview.github.io/?https://github.com/mrc-ide/global-lmic-reports-orderly/blob/master/vignettes/dide_run.html). 

For any troubleshooting help, please see the `FAQ.Rmd` vignette in the vignettes folder 
(also available to viewe [here](https://htmlpreview.github.io/?https://github.com/mrc-ide/global-lmic-reports-orderly/blob/master/vignettes/faq.html))
