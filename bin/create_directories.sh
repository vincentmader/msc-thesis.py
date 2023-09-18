#!/bin/sh

BLUE="\033[0;34m"
NC="\033[0m"

directories=(
    "../lib"
    "../out"
)

msg="Creating directories..."
echo "\n${BLUE}${msg}${NC}"
for dir in "${directories[@]}"; do
    echo "\`${dir}\`"
    if [ -d "${dir}" ]; then
        printf " -> Already exists."
    else
        mkdir "${dir}"
    fi
done
