import pandas as pd


def infer_and_add_platemap(metadata_df):
    plates_to_remove = (
        (
            metadata_df.groupby("Metadata_Plate")
            .Metadata_JCP2022.apply(lambda x: x.nunique())
            .reset_index()
        )
        .query("Metadata_JCP2022 == 1")["Metadata_Plate"]
        .tolist()
    )

    metadata_df = metadata_df.query("Metadata_Plate not in @plates_to_remove")

    plates = metadata_df.Metadata_Plate.unique().tolist()
    platemap_df = pd.DataFrame({"Metadata_Plate": plates, "Metadata_Platemap": None})
    id = 1
    for i in range(len(plates) - 1):
        plate_1_name = plates[i]
        if (
            platemap_df.query("Metadata_Plate == @plate_1_name")["Metadata_Platemap"]
            .isnull()
            .values[0]
        ):
            platemap_name = f"platemap_{id}"
            id += 1
            platemap_df.loc[
                platemap_df.Metadata_Plate == plate_1_name, "Metadata_Platemap"
            ] = platemap_name
            plate_1_df = metadata_df.query("Metadata_Plate == @plate_1_name")[
                ["Metadata_Well", "Metadata_Symbol"]
            ]
            for j in range(i + 1, len(plates)):
                plate_2_name = plates[j]
                if (
                    platemap_df.query("Metadata_Plate == @plate_2_name")[
                        "Metadata_Platemap"
                    ]
                    .isnull()
                    .values[0]
                ):
                    plate_2_df = metadata_df.query("Metadata_Plate == @plate_2_name")[
                        ["Metadata_Well", "Metadata_Symbol"]
                    ]
                    merged_plate = pd.merge(
                        plate_1_df,
                        plate_2_df,
                        on=["Metadata_Well", "Metadata_Symbol"],
                        how="inner",
                    )
                    if len(merged_plate) > 0.95 * len(plate_2_df):
                        platemap_df.loc[
                            platemap_df.Metadata_Plate == plate_2_name,
                            "Metadata_Platemap",
                        ] = platemap_name

    metadata_df = metadata_df.merge(platemap_df, on="Metadata_Plate", how="left")

    return metadata_df
