#!/bin/bash

find . -type f -name "*.txt" | while read -r file; do
    dir=$(dirname "$file")
    base=$(basename "$file")

    newbase="${base//scalevar_01/scalevar_1}"

    if [[ "$base" != "$newbase" ]]; then
        echo "Renaming '$file' -> '$dir/$newbase'"
        mv -- "$file" "$dir/$newbase"
    fi
done

echo "Done."
