import sys
import argparse

parser = argparse.ArgumentParser(
    description="Use the spread and the over/under to compute expected scores of two teams."
)
parser.add_argument(
    "combined",
    help="The expected combined score or over/under (O/U).",
    type=float,
)
parser.add_argument(
    "spread",
    help="The expected spread in scores.",
    type=float,
)
args = parser.parse_args()

ou = args.combined
spread = args.spread

print(f"Team 1 score: {spread + (ou - spread) / 2}")
print(f"Team 2 score: {(ou - spread) / 2}")
