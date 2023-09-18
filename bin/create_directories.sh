#!/bin/sh

BLUE="\033[0;34m"
NC="\033[0m"

directories=(
    "../out"
)

msg="Creating directories..."
echo "\n${BLUE}${msg}${NC}"
for dir in "${directories[@]}"; do
    printf " - ${dir}"
    if [ -d "${dir}" ]; then
        echo " -> Already exists."
    else
        mkdir "${dir}"
        echo
    fi
done
