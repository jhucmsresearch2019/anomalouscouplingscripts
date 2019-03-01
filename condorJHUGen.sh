#!/bin/bash

set -euo pipefail

cd $CMSSW_BASE
eval $(scram ru -sh) #this is equivalent to cmsenv
cd your JHUGen folder

mkdir -p outputfolder
./JHUGen (arguments) DataFile=outputfolder/$1.lhe Seed=123456$1
