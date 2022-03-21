#!/usr/bin/env bash

set -e
script_dir=`dirname "$0"`
cd "${script_dir}"/../
./scripts/install-r-dependencies.R

export PATH=${PATH}:`pwd`/

# No obvious other dependencies for now
