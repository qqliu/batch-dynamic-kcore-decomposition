#!/bin/bash

command -v make >/dev/null 2>&1 || { sudo apt install make; }
sudo apt-get update
command -v g++ >/dev/null 2>&1 || { sudo apt install g++; }

git submodule update --init --recursive
echo "Setup done!"
