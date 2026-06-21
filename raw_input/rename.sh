#!/bin/bash

find . -depth -name '*nNNPDF30*' | while read -r file; do
    newfile="${file//nNNPDF30/NNPDF30}"
    mv -- "$file" "$newfile"
done
