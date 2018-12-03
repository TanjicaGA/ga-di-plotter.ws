#!/bin/sh

WATCHED_DIR=${1-./}
PORT=${2-5000}

start () {
  R -e "shiny::runApp('$WATCHED_DIR', port = $PORT, host='0.0.0.0')" &
  PID=$!
}

start

inotifywait -mr $WATCHED_DIR --format '%e %f' \
  -e modify -e delete -e move -e create \
  | while read event file; do

  echo $event $file

  kill $PID
  start

done
