#!/bin/bash
#
# setup and configure WhatsHap
#
base="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd )"

whatshap=$base/whatshap/venv/bin/whatshap

# WhatsHap
#----------------------------------------------------------------------
echo "setting up WhatsHap .."

# clone git repo
function clone_whatshap {
    cd $base
    echo "cloning WhatsHap git repository .."
    git clone --branch matrix https://bitbucket.org/whatshap/whatshap.git
}

# setup virtual environment
function venv_whatshap {
    cd $base
    echo "setting up virtual environment .."
    cd whatshap
    virtualenv -p python3 venv
    venv/bin/pip3 install Cython nose tox
    venv/bin/pip3 install -e .
}

if [ -e $whatshap ] ; then
    echo "WhatsHap executable already exists"
else
    clone_whatshap
    venv_whatshap
fi

echo "setup of WhatsHap is complete."

# finish
#----------------------------------------------------------------------

echo
echo "setup of WhatsHap has been successful." 
echo
echo "WhatsHap can now be executed as:"
echo
echo "   $whatshap"
echo
