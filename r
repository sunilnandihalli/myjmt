#!/bin/bash -v
mkdir -p build && cd build && cmake ../ && make && ./myjmt
