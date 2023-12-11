#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: $0 <project-name>"
    exit 1
fi

project=$1

mkdir -p .build

pdflatex -output-directory .build -jobname $project $project.tex

mv .build/$project.pdf .
