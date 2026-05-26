find . -type f -name '*cbarW*' | while read -r file; do
    newfile=$(echo "$file" | sed 's/cbarW/cW/g')
    
    if [ "$file" != "$newfile" ]; then
        echo "Renaming: $file -> $newfile"
        mv "$file" "$newfile"
    fi
done
