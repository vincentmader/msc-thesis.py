#!/bin/sh

BLUE="\033[0;34m"
NC="\033[0m"

msg="Cloning custom themes for Python library \`matplotlib\`..."
echo "\n${BLUE}${msg}${NC}"

path="../lib/mpl-styles"
if [ ! -d $path ]; then
    url="https://github.com/vincentmader/mpl-styles"
    git clone $url $path
else
    cd $path && git pull
fi
