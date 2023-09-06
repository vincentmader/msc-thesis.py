#!/bin/sh

url="https://github.com/vincentmader/msc-thesis.py"
path="./msc-thesis.py"

git clone "$url" "$path"

cd "$path" && make setup
