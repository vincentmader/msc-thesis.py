#!/bin/sh

clear

BLUE="\033[0;34m"
NC="\033[0m"

msg="Cloning repository from \"https://github.com/vincentmader/msc-thesis.py\"..."
echo "\n${BLUE}${msg}${NC}"

url="https://github.com/vincentmader/msc-thesis.py"
path="./msc-thesis.py"

git clone "$url" "$path"

cd "$path" && make setup
