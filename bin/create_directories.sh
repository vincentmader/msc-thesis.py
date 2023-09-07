#!/bin/sh

BLUE="\033[0;34m"
NC="\033[0m"

directories=(
    "../lib"
    "../out"
)

for dir in "${directories[@]}"; do
    msg="Creating directory \`${dir}\`..."
    echo "\n${BLUE}${msg}${NC}"
    if [ -d "${dir}" ]; then
        echo "Already exists."
    else
        mkdir "${dir}"
    fi
done
