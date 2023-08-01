#!/bin/bash

input_file=$1
output_folder=$2

mkdir -p "${output_folder}"
mkdir -p "${output_folder}/logs"

fastqc -t 4 -o "${output_folder}" "${input_file}" > "${output_folder}/logs/$(basename "${input_file}").log" 2>&1
