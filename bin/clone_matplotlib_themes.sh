#!/bin/sh

path="../lib/mpl-styles"
if [ ! -d $path ]; then
    url="https://github.com/vincentmader/mpl-styles"
    git clone $url $path
else
    cd $path && git pull
fi
