#!/bin/sh

path="../lib/coag_py"
if [ ! -d $path ]; then
    # url="git@github.com:vincentmader/coag_py.git"
    url="https://github.com/vincentmader/coag_py"
    git clone $url $path
else
    cd $path && git pull
fi
