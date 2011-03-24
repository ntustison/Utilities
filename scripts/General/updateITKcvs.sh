# request Bourne shell as shell for job
#$ -S /bin/sh

cd /Users/nick/pkg/InsightToolkit/src/InsightToolkitCVS/Insight/
cvs update -d
cd /Users/nick/pkg/InsightToolkit/bin/
cmake /Users/nick/pkg/InsightToolkit/src/InsightToolkitCVS/Insight/
make

# chmod -R ugo+rwx /Users/nick/pkg/InsightToolkit/src/InsightToolkitCVS/InsightCopy/
# /bin/rm -r /Users/nick/pkg/InsightToolkit/src/InsightToolkitCVS/InsightCopy/
# cp -R /Users/nick/pkg/InsightToolkit/src/InsightToolkitCVS/Insight/ /Users/nick/pkg/InsightToolkit/src/InsightToolkitCVS/InsightCopy/
# cd /Users/nick/pkg/InsightToolkit/src/InsightToolkitCVS/InsightCopy/
# find ./ -name "CVS" | xargs rm -r
# find ./ -name ".cvsignore" | xargs rm -r
# cvs import -I ! -m "Nick's ITK" Insight Nick start
