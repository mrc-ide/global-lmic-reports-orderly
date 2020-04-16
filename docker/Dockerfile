FROM vimc/orderly:master

RUN apt-get update && apt-get -y install \
        libv8-dev

COPY bin /usr/local/bin/

RUN install_packages --repo=https://ncov-ic.github.io/drat \
        dde \
        dplyr \
        ggplot2 \
        ggrepel \
        knitr \
        lubridate \
        magrittr \
        odin \
        purrr \
        readxl \
        tidyr

RUN install_packages --repo=https://ncov-ic.github.io/drat \
        odin.js \
        squire
