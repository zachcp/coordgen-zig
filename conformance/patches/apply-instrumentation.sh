#!/bin/sh
set -eu

source_file=$1
patch_file=$2
output_file=$3

patch --batch --silent --output="$output_file" "$source_file" < "$patch_file"
