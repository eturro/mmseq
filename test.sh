#!/bin/bash

set -ex

pushd src
make clean
make
popd
