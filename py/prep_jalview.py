#!/usr/bin/env python3
"""Remap Jalview features from full-length sequence coords onto a prepped MSA.

Reads the original mixed-case MSA and the msa_prep map (muscle -map /
masm_train default <output>.map). A residue at orig column C is in the
seed iff C appears in the map; the ungapped seed index is this sequence's
letter count in mapped columns before C.

Prep ungapped coordinates are 1-based indices into the ungapped seed
sequence, including letters in kept insert/spacer columns, so Jalview
colors line up. Residues in discarded orig columns are omitted.

Supports classic Jalview feature lines:
  description  seq_id  group  start  end  featureType

and Muscle-style lines:
  -  seq_id  seq_index  start  end  featureType

Usage:
  python py/prep_jalview.py MOTIFS.jalview ORIG.afa ORIG.map -o OUT.jalview
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Optional, Tuple


def trunc_label(label: str) -> str:
	"""Match muscle default trunc_label: cut at first whitespace, '|', or '/'."""
	n = len(label)
	for i, c in enumerate(label):
		if c in " \t\r\n|/":
			n = i
			break
	return label[:n]


def label_key(label: str) -> str:
	return trunc_label(label).upper()


def is_gap(c: str) -> bool:
	return c in "-."


def read_fasta(path: str) -> Tuple[Dict[str, str], Dict[str, str]]:
	seqs: Dict[str, str] = {}
	canonical: Dict[str, str] = {}
	label: Optional[str] = None
	parts: List[str] = []

	def flush() -> None:
		nonlocal label, parts
		if label is None:
			return
		key = label_key(label)
		seq = "".join(parts)
		if key in seqs:
			raise SystemExit("%s: duplicate FASTA label %s" % (path, label))
		seqs[key] = seq
		canonical[key] = trunc_label(label)
		parts = []

	for line in open(path, encoding="utf-8", errors="replace"):
		line = line.rstrip("\n\r")
		if line.startswith(">"):
			flush()
			label = line[1:]
			parts = []
			continue
		if label is None:
			continue
		parts.append(line)
	flush()
	if len(seqs) == 0:
		raise SystemExit("%s: empty FASTA" % path)
	return seqs, canonical


def parse_region_fields(fields: List[str]) -> Tuple[str, Dict[str, str]]:
	kind = fields[0]
	i = 1
	if i < len(fields) and fields[i].isdigit():
		nxt = fields[i + 1] if i + 1 < len(fields) else ""
		if nxt in ("orig_lo", "orig_hi", "seed_lo", "seed_hi", "orig_cols"):
			i += 1
	d: Dict[str, str] = {}
	while i + 1 < len(fields):
		d[fields[i]] = fields[i + 1]
		i += 2
	return kind, d


def parse_msa_prep_map(path: str) -> Tuple[int, List[int]]:
	kept: List[int] = []
	orig_cols: Optional[int] = None
	seed_cols: Optional[int] = None
	saw_header = False
	for line in open(path, encoding="utf-8", errors="replace"):
		line = line.rstrip("\n\r")
		if not line:
			continue
		fields = line.split("\t")
		if fields[0] == "msa_prep_map":
			if len(fields) < 2 or fields[1] != "1":
				raise SystemExit("%s: unsupported map version %r" % (path, line))
			saw_header = True
			continue
		if not saw_header:
			raise SystemExit("%s: expected msa_prep_map header, got %r" % (path, line))
		if fields[0] == "orig_cols":
			orig_cols = int(fields[1])
			continue
		if fields[0] == "seed_cols":
			seed_cols = int(fields[1])
			continue
		kind, d = parse_region_fields(fields)
		if kind == "discard":
			continue
		if kind == "block":
			lo = int(d["orig_lo"])
			hi = int(d["orig_hi"])
			if hi < lo:
				raise SystemExit("%s: block orig_hi < orig_lo" % path)
			kept.extend(range(lo, hi + 1))
			continue
		if kind == "spacer":
			cols = d.get("orig_cols", "")
			if cols == "":
				raise SystemExit("%s: spacer missing orig_cols" % path)
			kept.extend(int(x) for x in cols.split(",") if x != "")
			continue
		raise SystemExit("%s: unknown map record %r" % (path, fields[0]))

	if not saw_header:
		raise SystemExit("%s: missing msa_prep_map header" % path)
	if orig_cols is None or seed_cols is None:
		raise SystemExit("%s: missing orig_cols/seed_cols" % path)
	if len(kept) != seed_cols:
		raise SystemExit(
			"%s: listed %d orig columns, seed_cols=%d" %
			(path, len(kept), seed_cols)
		)
	return orig_cols, kept


def map_one_sequence(row: str, kept_orig_cols: List[int]) -> Dict[int, int]:
	kept_set = set(kept_orig_cols)
	pos_map: Dict[int, int] = {}
	orig_pos = 0
	seed_pos = 0
	for col, c in enumerate(row):
		if is_gap(c):
			continue
		if col in kept_set:
			pos_map[orig_pos] = seed_pos
			seed_pos += 1
		orig_pos += 1
	return pos_map


def load_orig_map(
	msa_path: str, map_path: str
) -> Tuple[Dict[str, Dict[int, int]], Dict[str, str]]:
	seqs, canonical = read_fasta(msa_path)
	orig_cols, kept = parse_msa_prep_map(map_path)
	widths = {len(s) for s in seqs.values()}
	if len(widths) != 1:
		raise SystemExit("%s: uneven sequence lengths" % msa_path)
	width = widths.pop()
	if width != orig_cols:
		raise SystemExit(
			"%s: MSA width %d != map orig_cols %d" %
			(msa_path, width, orig_cols)
		)

	out: Dict[str, Dict[int, int]] = {}
	for key, row in seqs.items():
		out[key] = map_one_sequence(row, kept)
	return out, canonical


def contiguous_runs(positions: List[int]) -> List[Tuple[int, int]]:
	if not positions:
		return []
	positions = sorted(set(positions))
	runs: List[Tuple[int, int]] = []
	lo = hi = positions[0]
	for p in positions[1:]:
		if p == hi + 1:
			hi = p
		else:
			runs.append((lo, hi))
			lo = hi = p
	runs.append((lo, hi))
	return runs


def map_range(
	pos_map: Dict[int, int], start: int, end: int
) -> List[Tuple[int, int]]:
	if start > end:
		start, end = end, start
	mapped: List[int] = []
	for p in range(start, end + 1):
		orig0 = p - 1
		if orig0 in pos_map:
			mapped.append(pos_map[orig0])
	runs0 = contiguous_runs(mapped)
	return [(lo + 1, hi + 1) for lo, hi in runs0]


def parse_feature_fields(
	fields: List[str],
) -> Optional[Tuple[str, str, str, int, int, str]]:
	if len(fields) < 6:
		return None
	try:
		start = int(fields[3])
		end = int(fields[4])
	except ValueError:
		return None
	return (fields[0], fields[1], fields[2], start, end, fields[5])


def remap_jalview(
	jalview_path: str,
	prep: Dict[str, Dict[int, int]],
	canonical: Dict[str, str],
	out_path: str,
) -> None:
	missing_labels: set = set()
	skipped_empty = 0
	written_features = 0

	with open(jalview_path, "r", encoding="utf-8", errors="replace") as fin, open(
		out_path, "w", encoding="utf-8", newline="\n"
	) as fout:
		for line_num, line in enumerate(fin, 1):
			raw = line.rstrip("\n\r")
			if not raw:
				fout.write("\n")
				continue

			upper = raw.upper()
			if upper.startswith("STARTGROUP") or upper.startswith("ENDGROUP"):
				fout.write(raw + "\n")
				continue

			fields = raw.split("\t")
			if len(fields) == 1:
				fields = raw.split()

			if len(fields) == 2 and fields[0] != "-":
				fout.write(raw + "\n")
				continue

			parsed = parse_feature_fields(fields)
			if parsed is None:
				fout.write(raw + "\n")
				continue

			desc, seq_id, group, start, end, feat_type = parsed
			key = label_key(seq_id)
			if key not in prep:
				if key not in missing_labels:
					print(
						"warning: label %r not in orig MSA; skipping features" %
						seq_id,
						file=sys.stderr,
					)
					missing_labels.add(key)
				continue

			runs = map_range(prep[key], start, end)
			if not runs:
				skipped_empty += 1
				continue

			out_label = canonical.get(key, trunc_label(seq_id))
			for lo, hi in runs:
				fout.write(
					"%s\t%s\t%s\t%d\t%d\t%s\n" %
					(desc, out_label, group, lo, hi, feat_type)
				)
				written_features += 1

	msg = "wrote %d feature range(s) to %s" % (written_features, out_path)
	if skipped_empty:
		msg += " (%d unmapped ranges skipped)" % skipped_empty
	print(msg, file=sys.stderr)


def main(argv: Optional[List[str]] = None) -> int:
	ap = argparse.ArgumentParser(
		description="Map Jalview features onto prepped seed-MSA coordinates"
	)
	ap.add_argument("jalview", help="input Jalview features (full-length coords)")
	ap.add_argument("msa", help="original mixed-case MSA that was prepped")
	ap.add_argument("map", help="msa_prep_map TSV from muscle -map")
	ap.add_argument("-o", "--output", required=True, help="output Jalview features")
	args = ap.parse_args(argv)

	if not os.path.isfile(args.msa):
		raise SystemExit("orig MSA not found: %s" % args.msa)
	if not os.path.isfile(args.map):
		raise SystemExit("map file not found: %s" % args.map)

	prep, canonical = load_orig_map(args.msa, args.map)
	if not prep:
		raise SystemExit("no sequences in orig MSA %s" % args.msa)
	remap_jalview(args.jalview, prep, canonical, args.output)
	return 0


if __name__ == "__main__":
	sys.exit(main())
