#!/bin/sh

path="../src/solver/kees_solvers"
if [ ! -d $path ]; then
    url="git@github.com:vincentmader/coag_py.git"
    git clone $url $path
else
    cd $path && git pull
fi
