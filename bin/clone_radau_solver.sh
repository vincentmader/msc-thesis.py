#!/bin/sh

path="../src/solver/kees_solvers"
if [ ! -d $path ]; then
    url="https://github.com/vincentmader/coag_py"
    git clone $url $path
else
    cd $path && git pull
fi
