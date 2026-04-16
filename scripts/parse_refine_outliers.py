"""Parse augur refine output to extract strains with dates reset as outliers."""

import re
import sys

import pandas as pd

SENTINEL = "Inferred a time resolved phylogeny using TreeTime:"
OUTLIER_HEADER = "the following tips have been marked as outliers"


def format_date(s):
    """Format a date string (single float or tuple of floats) to 2 decimal places."""
    s = s.strip()
    if s.startswith("(") and s.endswith(")"):
        parts = s[1:-1].split(",")
        return "(" + ", ".join(f"{float(p):.2f}" for p in parts) + ")"
    return f"{float(s):.2f}"


def main():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")  # noqa: F821

    bad_dates_action = snakemake.params.bad_dates_action  # noqa: F821
    if bad_dates_action not in ("ignore", "error"):
        raise ValueError(f"Invalid {bad_dates_action=}, must be 'ignore' or 'error'")

    with open(snakemake.input.refine_output) as f:  # noqa: F821
        text = f.read()

    # Validate augur refine completed successfully
    if SENTINEL not in text:
        raise ValueError(
            f"Could not find expected sentinel '{SENTINEL}' in augur refine output. "
            "The output format may have changed or the run may be incomplete."
        )

    outlier_start = text.find(OUTLIER_HEADER)
    sentinel_pos = text.find(SENTINEL)

    entries = []
    if outlier_start >= 0:
        block = text[outlier_start:sentinel_pos]
        header_end_marker = "have been reset:"
        header_end = block.find(header_end_marker)
        if header_end < 0:
            raise ValueError(
                "Found outlier header but not expected 'have been reset:' ending. "
                "The augur refine output format may have changed."
            )
        block = block[header_end + len(header_end_marker) :]

        # Collapse continuation lines (start with whitespace) onto previous line
        lines = []
        for line in block.split("\n"):
            if not line.strip():
                continue
            if re.match(r"\d", line):
                lines.append(line)
            elif lines:
                lines[-1] += " " + line.strip()

        # Parse: <timestamp>\t<accession>, input date: <date>, apparent date: <date>
        entry_re = re.compile(r"(\S+), input date: (.+?), apparent date: (\S+)")
        for line in lines:
            match = entry_re.search(line)
            if match:
                entries.append(
                    {
                        "accession": match.group(1),
                        "metadata_date": format_date(match.group(2)),
                        "inferred_date": format_date(match.group(3)),
                    }
                )

        if not entries:
            raise ValueError(
                "Found outlier header but could not parse any entries. "
                "The augur refine output format may have changed."
            )

    # Safety: "outlier" should not appear outside the known block
    if outlier_start >= 0:
        text_outside = text[:outlier_start] + text[sentinel_pos:]
    else:
        text_outside = text
    if "outlier" in text_outside.lower():
        raise ValueError(
            "Found unexpected mention of 'outlier' outside the known outlier block "
            "in augur refine output. The output format may have changed."
        )

    if entries and bad_dates_action == "error":
        msg = (
            "The following accessions in keep-ids had their dates reset as outliers "
            "by augur refine, and bad_dates_in_keep_accessions_action is 'error':\n"
        )
        for e in entries:
            msg += (
                f"  {e['accession']}: metadata_date={e['metadata_date']}, "
                f"inferred_date={e['inferred_date']}\n"
            )
        raise ValueError(msg)

    df = pd.DataFrame(entries, columns=["accession", "metadata_date", "inferred_date"])
    df.to_csv(snakemake.output.tsv, sep="\t", index=False)  # noqa: F821
    print(f"Wrote {len(df)} outlier entries to {snakemake.output.tsv}")  # noqa: F821


if __name__ == "__main__":
    main()
