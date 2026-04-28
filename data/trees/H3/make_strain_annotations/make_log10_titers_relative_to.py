import marimo

__generated_with = "0.20.4"
app = marimo.App(width="full")


@app.cell
def _():
    import os

    import numpy

    import pandas as pd

    input_csv = "virus_titer_aggregated.csv"
    output_csv = "log10_titers_relative_to.csv"

    norm_to = ["293-noSA", "293-SA23", "293-SA26"]

    print(f"Reading input from {input_csv=}")
    input_titers = pd.read_csv(input_csv)
    print(f"Found the following columns:\n{input_titers.columns.tolist()}")
    assert len(input_titers) == len(input_titers.groupby(["strain", "cell_line"]))

    subtype = os.path.basename(os.path.abspath("../"))
    print(f"Filtering for {subtype=}")
    n_pre = len(input_titers)
    input_titers = input_titers[input_titers["subtype"] == subtype]
    print(f"Only keeping the {len(input_titers)} of {n_pre} rows that have {subtype=}")

    strains_to_drop = {"A/swine/Changhua/199-3/2000", "A/WSN/1933-H141Y"}
    for strain_to_drop in strains_to_drop:
        if strain_to_drop in set(input_titers["strain"]):
            n_pre = len(input_titers)
            input_titers = input_titers[input_titers["strain"] != strain_to_drop]
            print(f"Dropping rows for {strain_to_drop=}: kept {len(input_titers)} of {n_pre} rows")

    print("Converting to wide form log10 titers")
    titers_wide = input_titers.pivot_table(
        index="strain",
        columns="cell_line",
        values="mean_RLUperuL",
    )
    assert titers_wide.notnull().all().all()
    assert (titers_wide > 0).all().all()
    log10_titers_wide = numpy.log10(titers_wide)

    print(f"Now computing log10 titers relative to {norm_to}")
    assert set(titers_wide.columns).issuperset(norm_to)
    normed_data = {}
    for norm_col in norm_to:
        for col in titers_wide.columns:
            if col != norm_col:
                normed_data[f"log10_{col}_over_{norm_col}_titer"] = log10_titers_wide[col] - log10_titers_wide[norm_col]
    log10_titers_relative_to = pd.DataFrame(normed_data)
    print(f"Computed the following columns:\n{log10_titers_relative_to.columns.tolist()}")
    log10_titers_relative_to.insert(0, "has_titers", "yes")
    print(f"Writing output to {output_csv}")
    log10_titers_relative_to.to_csv(output_csv, float_format="%.2f")
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
