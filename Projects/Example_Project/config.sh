#!/bin/bash
# Defines the paths a project's run scripts need. Source this from other
# scripts in this folder; don't run it directly.

set -e

SH_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# This project's own directory (config.sh lives directly in it here, not
# in a scripts/ subfolder)
PROJECT_DIR="$SH_DIR"

# The shared solver repo -- one build, many project directories point at it.
LS_IBM_CPP_DIR="/home/groups/ibattiat/sxia/LS_IBM/LS_IBM_sxia/LS_IBM_cpp"

# This project's config/input/output
CONFIG_DIR="$PROJECT_DIR/config"
INPUT_DIR="$PROJECT_DIR/input"
OUTPUT_DIR="$PROJECT_DIR/output"
