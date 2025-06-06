## Did you miss the previous step? The one about setting up Git and GitHub
rstudioapi::navigateToFile(usethis::proj_path("dev", "02_git_github_setup.R"))

## ***********************************************************
## Setup the core files for your Bioconductor-friendly package
## ***********************************************************

## Edit your package DESCRIPTION file
## Check http://r-pkgs.had.co.nz/description.html for details
## as well as http://contributions.bioconductor.org/description.html

## Check https://github.com/lcolladotor/biocthis/blob/devel/DESCRIPTION
## for an example.

## You'll at least want to edit the version to 0.99.0 (or lower) since that's
## the version number you will need to use with Bioconductor.

## You will also want to add the biocViews field, for example:
## biocViews: Software

## Many Bioconductor packages use the following license:
## license: Artistic-2.0

## You might want to add the Date field as well, which is used when creating
## the package citation information. Use the YYYY-MM-DD format. For example:
## Date: 2020-04-29

## This function sets all these defaults for you
biocthis::use_bioc_description()
## However, you still need to edit parts of it manually
rstudioapi::navigateToFile(usethis::proj_path("DESCRIPTION"))

## Create your README.Rmd file
biocthis::use_bioc_readme_rmd()
devtools::build_readme()

## Edit accordingly. You might want to install your package also using
## devtools::build() or the RStudio keyboard shortcut:
## macoS: command + shift + B
## Windows: control + shift + B

## Click on the `knit` button on your README.Rmd file to create the README.md
## file.

## Add a NEWS.md file
## See http://bioconductor.org/developers/package-guidelines/#news for more
## details about Bioconductor news files.
biocthis::use_bioc_news_md()

## Add information for users and contributors
biocthis::use_bioc_coc()
usethis::use_tidy_contributing()
biocthis::use_bioc_support()
biocthis::use_bioc_issue_template()
biocthis::use_bioc_feature_request_template()
biocthis::use_bioc_citation()

## Add badges to the README.Rmd file
usethis::use_lifecycle_badge("Experimental")
biocthis::use_bioc_badges()
## NOTE: If your Bioconductor package is an experiment, annotation or workflow
## package, you will need to edit the resulting badge!
usethis::use_github_actions_badge("check-bioc.yml")

## Enable using tests
usethis::use_testthat()
usethis::use_test("example_test") ## You need at least one test to run covr
usethis::use_coverage()

## Re-knit your README.Rmd file to update your README.md file
devtools::build_readme()

## Add a vignette template
pkg <- basename(usethis::proj_get())
biocthis::use_bioc_vignette(pkg, paste("Introduction to", pkg))

## Add a Bioconductor-friendly GitHub actions workflow to check your package
biocthis::use_bioc_github_action()
## You will need to go to https://github.com/eturkes/genefunnel/settings/actions
## to and enable:
## > Workflow permissions > Read and write permissions
## Then click save before you can continue with the instructions in this script.

## While biocthis::use_bioc_github_action() has arguments for all these options,
## if you missed them, and:
## * your package doesn't have testthat tests, use: has_testthat = 'false'
## * you don't want to run the covr step, use: run_covr = 'false'
## * you don't want to use pkgdown, use: run_pkgdown = 'false'
## * you don't want to build a docker image, use: run_docker = 'false'
rstudioapi::navigateToFile(usethis::proj_path(".github", "workflows", "check-bioc.yml"))

## Setup up your global git config
usethis::edit_git_config()
## Use the information that matches your GitHub account
## Example contents (4 space indentation):
# [user]
#     name = Your Full Name
#     email = your.email@somewhere.com
#

## Setup your pkgdown css file
biocthis::use_bioc_pkgdown_css()

## Install your package locally before building the documentation website
devtools::install()

## ************************* WARNING START *********************************
## WARNING: git commit before running this next function!
## Otherwise you can lose your work!!!
## ************************* WARNING END ***********************************
##
## Deploy with pkgdown at least once locally such that the automatic updates
## from GitHub actions will work. This creates the gh-pages branch in your
## GitHub repository in such a way that pkgdown will recognize it and be
## able to use it later.
pkgdown::deploy_to_branch() ## Check the WARNING above before running this!

## For more details about package development, check:
## * http://contributions.bioconductor.org/
## * https://r-pkgs.org/

## Move to the next step: updating your package code before a "git commit"
rstudioapi::navigateToFile(usethis::proj_path("dev", "04_update.R"))

## This template was made using https://lcolladotor.github.io/biocthis/
