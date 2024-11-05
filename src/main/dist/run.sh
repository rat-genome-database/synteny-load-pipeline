#!/usr/bin/env bash
# shell script to run Synteny Load pipeline
. /etc/profile

APPNAME=SyntenyLoadPipeline
APPDIR=/home/rgddata/pipelines/$APPNAME
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

cd $APPDIR

java -Dspring.config=$APPDIR/../properties/default_db2.xml \
    -Dlog4j.configurationFile=file://$APPDIR/properties/log4j2.xml \
    -jar lib/${APPNAME}.jar "$@" > $APPDIR/run.log 2>&1

mailx -s "[$SERVER] Synteny Load pipeline ok" mtutaj@mcw.edu < $APPDIR/logs/summary.log
