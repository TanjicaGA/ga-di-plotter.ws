#!/bin/sh

target=gamap:/var/www/shiny/DI-Plotter
source="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd $source;

rsync -cuv \
      "$source/renv.lock" "$target/"
