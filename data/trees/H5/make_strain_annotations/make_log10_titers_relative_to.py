import marimo

__generated_with = "0.16.5"
app = marimo.App(width="full")


@app.cell
def _():
    import numpy

    import pandas as pd

    input_csv = "virus_titer_QC_pass.csv"
    output_csv = "log10_titers_relative_to.csv"

    norm_to = ["293-SA23", "293-noSA"]

    print(f"Reading input from {input_csv=}")
    input_titers = pd.read_csv(input_csv)
    print(f"Found the following columns:\n{input_titers.columns.tolist()}")
    assert len(input_titers) == len(input_titers.groupby(["strain", "cell_line"]))

    print("Converting to wide form log10 titers")
    titers_wide = input_titers.pivot_table(
        index=["strain", "clade"],
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
