#!/bin/sh

BLUE="\033[0;34m"
NC="\033[0m"

msg="Cloning repository forked from \"https://github.com/dullemond/coag_py\"..."
echo "\n${BLUE}${msg}${NC}"

path="../lib/coag_py"
if [ ! -d $path ]; then
    # url="git@github.com:vincentmader/coag_py.git"
    url="https://github.com/vincentmader/coag_py"
    git clone $url $path
else
    cd $path && git pull
fi
