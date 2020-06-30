#!/bin/sh

WATCHED_DIR=${1-./}
PORT=${2-5000}

start () {
  R --slave -e "shiny::runApp('$WATCHED_DIR', port = $PORT, host='0.0.0.0')"
  PID=$!
}

start
