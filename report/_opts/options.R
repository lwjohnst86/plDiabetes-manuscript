##'
##' Code preamble
##' =============
##'

##' Load relevant libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(magrittr)
library(tidyr)
library(knitr)
library(pander)
library(geepack) ## For gee analysis
library(stringr)
library(broom)
library(grid)
library(lsmeans)
library(MuMIn)
library(ggthemes)
library(captioner)
library(showtext)

##' Source personal functions, or install (or update) custom functions
##' in a package (to install from github, needs devtools)
if (!require(devtools)) install.packages('devtools')
if (!require(rstatsToolkit)) devtools::install_github('lwjohnst86/rstatsToolkit')
library(rstatsToolkit)

##' Custom options for pander
panderOptions('table.split.table', Inf)
panderOptions('table.style', 'rmarkdown')
panderOptions('table.alignment.default',
              function(df) ifelse(sapply(df, is.numeric), 'right', 'left'))

##' Knitr global options
knitr::opts_chunk$set(warning = FALSE, echo = FALSE,
                      fig.width = 8, dpi = 300,
                      message = FALSE, dev = c('png', 'pdf', 'postscript'),
                      fig.showtext = TRUE)

##' For the table and figure references
options(tabcap.prefix = "Table", suptabcap.prefix = "Supplemental Table",
        tabcap.sep = ":", tabcap.prefix.highlight = "**",
        figcap.prefix = "Figure", supfigcap.prefix = "Supplemental Fig.",
        figcap.sep = ":", figcap.prefix.highlight = "**")
