# request Bourne shell as shell for job
#$ -S /bin/sh

cd /Users/nick/pkg/PICSL/ANTS/bin/
# svn update
# cd /Users/nick/pkg/PICSL/ANTS/bin/
# cmake /Users/nick/pkg/PICSL/ANTS/Examples/
ctest -D NightlyStart
ctest -D NightlyUpdate
ctest -D NightlyConfigure
ctest -D NightlyBuild
ctest -D NightlySubmit
ctest -D NightlyTest
#ctest -D NightlyCoverage
ctest -D NightlySubmit
#ctest -D NightlyMemCheck
ctest -D NightlySubmit
