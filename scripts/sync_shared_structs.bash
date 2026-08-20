#!/usr/bin/env bash
# Sync shared STRUCTS sources from reseek into muscle/src.
# Canonical copy lives in reseek; do not hard-link.
set -euo pipefail

RESEEK_SRC="${RESEEK_SRC:-$(dirname "$0")/../../reseek/src}"
MUSCLE_SRC="$(cd "$(dirname "$0")/../src" && pwd)"
LIST="$MUSCLE_SRC/shared_structs.files"

if [[ ! -d "$RESEEK_SRC" ]]; then
	echo "RESEEK_SRC not found: $RESEEK_SRC" >&2
	echo "Set RESEEK_SRC to reseek/src" >&2
	exit 1
fi
if [[ ! -f "$LIST" ]]; then
	echo "Missing $LIST" >&2
	exit 1
fi

copied=0
while IFS= read -r line || [[ -n "$line" ]]; do
	[[ -z "$line" || "$line" =~ ^# ]] && continue
	src="$RESEEK_SRC/$line"
	dst="$MUSCLE_SRC/$line"
	if [[ ! -f "$src" ]]; then
		echo "SKIP missing in reseek: $line" >&2
		continue
	fi
	cp -f "$src" "$dst"
	copied=$((copied + 1))
	echo "copied $line"
done < "$LIST"

# Host-bundled alphabet blob (not in shared_structs.files list)
if [[ -f "$RESEEK_SRC/alpha_collect.cpp" ]]; then
	cp -f "$RESEEK_SRC/alpha_collect.cpp" "$MUSCLE_SRC/alpha_collect.cpp"
	echo "copied alpha_collect.cpp"
	copied=$((copied + 1))
fi

echo "Done. $copied files synced into $MUSCLE_SRC"
