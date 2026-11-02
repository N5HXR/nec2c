#!/usr/bin/env python3
"""Add deterministic analysis cards to an Antenna Turtle NEC geometry export."""

import argparse
import hashlib
from pathlib import Path


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--frequency-mhz", type=float, required=True)
    parser.add_argument("--theta-count", type=int, default=37)
    parser.add_argument("--phi-count", type=int, default=73)
    parser.add_argument("--theta-step", type=float, default=5.0)
    parser.add_argument("--phi-step", type=float, default=5.0)
    arguments = parser.parse_args()
    if arguments.theta_count < 1 or arguments.phi_count < 1:
        parser.error("pattern counts must be positive")
    lines = arguments.source.read_text(encoding="ascii").splitlines()
    if not lines or lines[-1].strip() != "EN":
        parser.error("source must end with an EN card")
    source_digest = sha256(arguments.source)
    analysis_cards = [
        f"FR 0 1 0 0 {arguments.frequency_mhz:.12g}",
        (f"RP 0 {arguments.theta_count} {arguments.phi_count} 1000 "
         f"0.0 0.0 {arguments.theta_step:.12g} {arguments.phi_step:.12g} 0.0 0.0"),
    ]
    provenance = f"CM PERFORMANCE FIXTURE SOURCE SHA256 {source_digest}"
    output = "\n".join([provenance, *lines[:-1], *analysis_cards, "EN", ""])
    arguments.destination.parent.mkdir(parents=True, exist_ok=True)
    arguments.destination.write_text(output, encoding="ascii")
    print(f"{arguments.destination}: {sha256(arguments.destination)}")


if __name__ == "__main__":
    main()
