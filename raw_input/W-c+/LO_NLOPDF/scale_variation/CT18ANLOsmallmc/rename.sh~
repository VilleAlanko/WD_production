#!/bin/bash

find . -type f -name "*.txt" | while IFS= read -r file; do
    dir=$(dirname "$file")
    base=$(basename "$file")

    newbase="${base//CT18ANLOmcvar/CT18ANLO}"

    if [[ "$base" != "$newbase" ]]; then
        mv -- "$file" "$dir/$newbase"
        echo "Renamed: $file -> $dir/$newbase"
    fi
done
