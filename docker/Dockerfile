FROM vimc/orderly:master

RUN apt-get update && apt-get -y install \
        libcurl4-openssl-dev \
        libv8-dev \
        zlib1g-dev

COPY bin /usr/local/bin/

RUN install_packages --repo=https://ncov-ic.github.io/drat \
        V8 \
        curl \
        dde \
        dplyr \
        ggforce \
        ggplot2 \
        ggrepel \
        knitr \
        lubridate \
        magrittr \
        odin \
        plotly \
        purrr \
        readxl \
        tidyr

## Might make a base image below this that we update out of sync to
## speed this up?
RUN install_packages --repo=https://ncov-ic.github.io/drat \
        odin.js \
        squire
