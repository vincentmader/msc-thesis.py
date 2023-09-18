#!/bin/sh

clear
BLUE="\033[0;34m"
NC="\033[0m"

msg="Cloning repository from \"https://github.com/vincentmader/msc-thesis.py\"..."
echo "\n${BLUE}${msg}${NC}"

url_1="git@github.com:vincentmader/msc-thesis.py.git"
url_2="https://github.com/vincentmader/msc-thesis.py"
path="~/Desktop/./msc-thesis.py"

clone() {
    git clone --recurse-submodules "$1" "$2"
}
clone "$url_1" "$path" || clone "$url_2" "$path"

cd "$path" && make setup
