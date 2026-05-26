#!/bin/bash

# Recursively remove "M55__" from all file names
find . -type f -name '*M55__*' | while read -r file; do
    dir=$(dirname "$file")
    base=$(basename "$file")
    newbase="${base//M55__/}"
    
    # Only rename if the name actually changes
    if [[ "$base" != "$newbase" ]]; then
        mv "$file" "$dir/$newbase"
        echo "Renamed: $file -> $dir/$newbase"
    fi
done
