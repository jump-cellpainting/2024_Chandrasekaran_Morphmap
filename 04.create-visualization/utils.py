import pandas as pd
import numpy as np

def infer_platemap(df):
    plates = df.Plate.unique().tolist()
    platemap_df = pd.DataFrame({"Plate": plates, "Platemap": np.nan})
    id = 1
    for i in range(len(plates) - 1):
        plate_1_name = plates[i]
        if (
            platemap_df.query("Plate == @plate_1_name")["Platemap"]
            .isnull()
            .values.any()
        ):
            platemap_name = f"platemap_{id}"
            id += 1
            platemap_df.loc[platemap_df.Plate == plate_1_name, "Platemap"] = (
                platemap_name
            )
            plate_1 = df.query("Plate == @plate_1_name")[["Metadata_Well", "Symbol"]]
            for j in range(i + 1, len(plates)):
                plate_2_name = plates[j]
                if (
                    platemap_df.query("Plate == @plate_2_name")["Platemap"]
                    .isnull()
                    .values.any()
                ):
                    plate_2 = df.query("Plate == @plate_2_name")[
                        ["Metadata_Well", "Symbol"]
                    ]
                    merged_plate = pd.merge(
                        plate_1, plate_2, on=["Metadata_Well", "Symbol"], how="inner"
                    )
                    if len(merged_plate) > 0.95 * len(plate_2):
                        platemap_df.loc[
                            platemap_df.Plate == plate_2_name, "Platemap"
                        ] = platemap_name

    return platemap_df


def identify_control_plate(df):
    test = (
        df.groupby("Plate")
        .agg({"Metadata_JCP2022": "nunique"})
    )

    print(test)