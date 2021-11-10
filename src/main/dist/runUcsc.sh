. /etc/profile

APPNAME=SyntenyLoadPipeline
APPDIR=/home/rgddata/pipelines/$APPNAME
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

cd $APPDIR

java -Dspring.config=$APPDIR/../properties/default_db.xml \
    -Dlog4j.configuration=file://$APPDIR/properties/log4j.properties \
    -jar lib/${APPNAME}.jar -ucsc "$@" > $APPDIR/ucsc.log 2>&1

mailx -s "[$SERVER] Synteny Load pipeline OK" mtutaj@mcw.edu < $APPDIR/logs/summary.log
