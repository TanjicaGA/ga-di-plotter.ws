#!/bin/sh

target=toolbox:apps/PlateAnalyzer
source=$GA_PATH/Code/PlateAnalyzer

rsync -lrv \
      --delete "$source/" "$target/" \
      --exclude deploy.sh \
      --exclude http_redirect.py \
      --exclude .Rprofile
