#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

for file in $SCRIPT_DIR/*.ui
do
  outFile="$(echo "$file" | cut -f 1 -d '.').py"
  pyuic5 $file -o $outFile
done
