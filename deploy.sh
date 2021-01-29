#!/bin/sh

target=gamap:/var/www/shiny/DI-Plotter
source=DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd $source;

rsync -lcruv \
      --delete \
      --exclude deploy.sh \
      --exclude http_redirect.py \
      --exclude "renv/library" \
      --exclude .git \
      --exclude .gitignore \
      "$source/" "$target/"
