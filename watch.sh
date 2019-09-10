#!/bin/sh

WATCHED_DIR=${1-./}
PORT=${2-5000}

start () {
  R -e "shiny::runApp('$WATCHED_DIR', port = $PORT, host='0.0.0.0')" &
  PID=$!
}

start

inotifywait -r -m $WATCHED_DIR --format '%e %f' ~/git/R-packages/ga.software.dd \
  -e modify -e delete -e move -e create \
  | while read event file; do

  echo $event $file

  if [ $file != '.#'* ]; then
      kill $PID
      start
  fi

done
