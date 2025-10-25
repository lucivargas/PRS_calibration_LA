#!/bin/bash

# Script to check if all the dependencies required are installed.

# Track whether any checks fail
MISSING=0

## Check command-line tools
check_command() {
    if ! command -v "$1" &> /dev/null; then
        echo " Error: '$1' is not installed or not in PATH."
        MISSING=1
    else
        echo " Found: $1"
    fi
}

## Check R packages
check_r_package() {
    Rscript -e "if (!requireNamespace(\"$1\", quietly = TRUE)) { quit(status = 1) }" &> /dev/null
    if [ $? -ne 0 ]; then
        echo " Error: R package '$1' is not installed."
        MISSING=1
    else
        echo " Found R package: $1"
    fi
}

# ---- Perform checks ----
check_command bcftools
check_command plink2
check_command R
check_command yq
check_command parallel
check_r_package tidyverse
check_r_package argparse
check_r_package stringr
check_r_package doParallel
check_r_package foreach
check_r_package boot

if [ "$MISSING" -ne 0 ]; then
    echo "❌ One or more required tools/packages are missing. Check details in log." 
    exit 1
else
    echo "✅ All required programs and packages are available."
    exit 0
fi


## End
