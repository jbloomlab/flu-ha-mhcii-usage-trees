#!/usr/bin/env python3
"""Generate auspice config with color scales for phenotype scores.

Reads node data JSON to compute min/max values for each phenotype score,
then generates an auspice config with viridis color scales.
"""

import argparse
import json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phenotypes", nargs="+", required=True)
    parser.add_argument("--continuous-scale", nargs="+", required=True)
    parser.add_argument("--node-data", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.node_data) as f:
        node_data = json.load(f)

    colorings = []
    for phenotype in args.phenotypes:
        for suffix in ["total effect", "max magnitude effect"]:
            key = f"{phenotype} {suffix}"
            values = [
                attrs[key]
                for attrs in node_data["nodes"].values()
                if key in attrs and isinstance(attrs[key], (int, float))
            ]
            if not values:
                raise ValueError(f"No numeric values found for {key!r} in node data")
            vmin, vmax = min(values), max(values)
            num_colors = len(args.continuous_scale)
            scale = [
                [vmin + (vmax - vmin) * i / (num_colors - 1), color]
                for i, color in enumerate(args.continuous_scale)
            ]
            colorings.append(
                {
                    "key": key,
                    "title": key,
                    "type": "continuous",
                    "scale": scale,
                }
            )

        colorings.append(
            {
                "key": f"{phenotype} max magnitude mutation",
                "title": f"{phenotype} max magnitude mutation",
                "type": "categorical",
            }
        )

    with open(args.output, "w") as f:
        json.dump({"colorings": colorings}, f, indent=2)


if __name__ == "__main__":
    main()
