#!/usr/bin/env bash
# run_kofamscan.sh
# Usage: ./run_kofamscan.sh path/to/file.faa

FAA_FILE="$1"

/Users/krzysztofmrozik/tools/kofam_scan/exec_annotation \
  --config /Users/krzysztofmrozik/tools/kofam_scan/config.yml \
  --profile /Users/krzysztofmrozik/kofam_db/profiles/ \
  --ko-list /Users/krzysztofmrozik/kofam_db/ko_list \
  "$FAA_FILE" > "${FAA_FILE%.faa}_kofamscan.txt"

