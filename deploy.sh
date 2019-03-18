#!/bin/sh

target=gamap:apps/PlateAnalyzer
source=~/git/ga-di-plotter

cd $source;

rsync -lruv \
      --delete \
      --exclude deploy.sh \
      --exclude http_redirect.py \
      --exclude "packrat/lib*" \
      --exclude .git \
      --exclude .gitignore \
      "$source/" "$target/"
