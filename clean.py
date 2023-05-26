import argparse
import os
import pickle

import numpy as np
import pandas as pd
import yaml
from rdkit import Chem
from rdkit.Chem import Descriptors


def sanitize_smiles(smiles):
    try:
        smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    except:
        print(f"Invalid SMILES: {smiles}")
        smiles = np.nan
    if smiles == "":
        print("Invalid SMILES: ''")
        smiles = np.nan
    return smiles


def creating_single_molecule_directory(output_dir, name, years):
    """
    Open the files and build a full dictionary

    input:
        name of the smiles directory

    output:
        dictionary with all the smiles molecules and patents
    """
    subject_smiles_dictionary_temp = {}
    subject_smiles_dictionary = {}
    for year in years:
        try:
            with open(
                os.path.join(
                    output_dir, name, "subject_smiles_dictionary_" + name + "_" + str(year)
                ),
                "rb",
            ) as fp:
                subject_smiles_dictionary_temp = pickle.load(fp)
            print(str(year), "directory is", len(subject_smiles_dictionary_temp), "long.")

            for k, v in subject_smiles_dictionary_temp.items():
                subject_smiles_dictionary.setdefault(k, set())
                subject_smiles_dictionary[k] |= set(v)
        except Exception as e:
            print("did not work", e)
            continue
    print("Length of", name, "directory is", len(subject_smiles_dictionary))
    return subject_smiles_dictionary


def cleaning_data_set(directory, charged_only=False, neutral_only=False):
    """
    Creates a pandas dataframe for easier overlooking of files.
    Also selects only the molecules which have a charge in the smiles
    representation of a molecule.

    input: smiles molecule python dictionary with patent names

    output: pandas DF with unique molecules and all the patents
        that have these molecules included in them
    """
    print("Size of directory before selecting charged molecules:", len(directory))
    smiles_DF = pd.DataFrame.from_dict(directory, orient="index")
    smiles_DF = smiles_DF.reset_index()
    smiles_DF = smiles_DF.rename(columns={"index": "smiles"})
    # select molecules with positive or negative charge
    smiles_DF["indexes_pos"] = smiles_DF.smiles.str.find("+")
    smiles_DF["indexes_neg"] = smiles_DF.smiles.str.find("-")
    if charged_only:
        # drop all the empty columns which are created when uncharged molecules are dropped
        df_charged = smiles_DF[(smiles_DF.indexes_pos >= 0) | (smiles_DF.indexes_neg >= 0)]
        df_charged = df_charged.dropna(axis=1, how="all")
        df_charged = df_charged.drop(["indexes_pos", "indexes_neg"], axis=1)
        print("Size of directory after selecting charged molecules:", df_charged.shape[0])
        return df_charged
    elif neutral_only:
        # drop all the empty columns which are created when charged molecules are dropped
        df_neutral = smiles_DF[(smiles_DF.indexes_pos == -1) & (smiles_DF.indexes_neg == -1)]
        df_neutral = df_neutral.dropna(axis=1, how="all")
        df_neutral = df_neutral.drop(["indexes_pos", "indexes_neg"], axis=1)
        print("Size of directory after selecting neutral molecules:", df_neutral.shape[0])
        return df_neutral


def cleaning_molecules_from_subst(df_molecules):
    """
    Cleaning the molecules from the different
    R (* in the dataset), Y substituents
    Changing R into ethylene group and Y into oxygen
    Also, dropping molecules which stayed invalid after change

    input:
        Pandas dataframe with smiles molecules

    output:
        Cleaned database with only chemically valid molecules
    """
    # Manually changing syntax of smiles molecules to ethylene groups using regex
    df_clean = df_molecules.copy()
    print('Substituting "*" and "Y" substituents with ethyl and O groups.')
    df_clean.smiles = df_clean.smiles.str.replace(r"\[\d+\*\]", "CC")
    df_clean.smiles = df_clean.smiles.str.replace(r"(\d\d\*\+)", "CC+")
    df_clean.smiles = df_clean.smiles.str.replace(r"(\d\*\+)", "CC+")
    df_clean.smiles = df_clean.smiles.str.replace(r"(\*\+)", "CC+")
    df_clean.smiles = df_clean.smiles.str.replace(r"\[\d+\*+:\d\]", "CC")
    df_clean.smiles = df_clean.smiles.str.replace(r"\[\*+:\d\]", "CC")
    df_clean.smiles = df_clean.smiles.str.replace(r"\*", "CC")
    df_clean.smiles = df_clean.smiles.str.replace("Y", "O")

    # Sanitizing molecules, checking for and dropping invalid molecules
    print("Sanitizing molecules, checking for invalid molecules...")
    df_clean["smiles"] = df_clean["smiles"].apply(lambda x: sanitize_smiles(x))
    print(f"Removing {df_clean.smiles.isna().sum()} invalid molecules...")
    df_clean.dropna(subset=["smiles"], inplace=True)
    print("Done.")

    # Collect all the patents into a single column
    df_clean["patents"] = df_clean[df_clean.columns[1:]].apply(
        lambda x: ",".join(x.dropna().astype(str)), axis=1
    )
    return df_clean


def select_molecules_by_mw(df_molecules, min_mw=0, max_mw=1e4):
    """
    Selects molecules by their molecular weight (MW) range.

    input:
        Pandas dataframe with smiles molecules
        min_mw: minimum molecular weight
        max_mw: maximum molecular weight

    output:
        Pandas dataframe with molecules within the MW range
    """
    print(f"Selecting molecules with MW between {min_mw} and {max_mw}...")
    df_molecules["MW"] = df_molecules["smiles"].apply(
        lambda x: Descriptors.MolWt(Chem.MolFromSmiles(x))
    )
    df_molecules = df_molecules[(df_molecules["MW"] >= min_mw) & (df_molecules["MW"] <= max_mw)]
    print(f"Done. {df_molecules.shape[0]} molecules selected.")
    return df_molecules


def get_parser():
    parser = argparse.ArgumentParser(description="Downloads patent files from USPTO website")
    parser.add_argument(
        "--years",
        type=str,
        nargs="+",
        required=True,
        help="Year(s) of patent files to download (or 'all')",
    )
    parser.add_argument(
        "--naming",
        type=str,
        required=True,
        help="Name of query (location where results will be stored)",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default=".",
        help="Parent directory of 'args.naming'",
    )
    parser.add_argument(
        "--charged_only",
        default=False,
        action="store_true",
        help="Whether to include only charged molecules in the dataset",
    )
    parser.add_argument(
        "--neutral_only",
        default=False,
        action="store_true",
        help="Whether to include only neutral molecules in the dataset",
    )
    parser.add_argument(
        "--mw_min",
        type=float,
        default=0,
        help="Minimum molecular weight to accept",
    )
    parser.add_argument(
        "--mw_max",
        type=float,
        default=1e4,
        help="Maximum molecular weight to accept",
    )
    return parser


def main(args):

    if args.years == ["all"]:
        print("Preparing to search chemistry patents from 2005 to 2023...")
        years = list(map(str, range(2005, 2023)))
        # TODO: previous step (search) doesn't currently work for 2001-2004 because of different XML structure
    else:
        print("Preparing to search chemistry patents from", ", ".join(args.years), "...")
        years = args.years

    subject_smiles_dictionary_abstract = creating_single_molecule_directory(
        args.output_dir, args.naming, years
    )

    smiles_DF = pd.DataFrame.from_dict(subject_smiles_dictionary_abstract, orient="index")
    smiles_DF = smiles_DF.reset_index()
    smiles_DF = smiles_DF.rename(columns={"index": "smiles"})

    if args.charged_only:
        df = cleaning_data_set(subject_smiles_dictionary_abstract, charged_only=True)
    elif args.neutral_only:
        df = cleaning_data_set(subject_smiles_dictionary_abstract, neutral_only=True)
    else:
        df = smiles_DF

    df_cleaned = cleaning_molecules_from_subst(df)

    if args.mw_min > 0 or args.mw_max < 1e4:
        df_cleaned = select_molecules_by_mw(df_cleaned, args.mw_min, args.mw_max)

    df_cleaned_smiles_patents = df_cleaned[["smiles", "patents"]].copy()

    # Create yaml file with SMILES and patents
    smiles_dict_filelist = []
    for ind in df_cleaned_smiles_patents.index:
        smiles_dict_filelist.append(
            {
                "smiles": df_cleaned_smiles_patents.smiles.loc[ind],
                "reference": [df_cleaned_smiles_patents.patents.loc[ind]],
            }
        )
    with open(
        os.path.join(
            args.output_dir,
            args.naming,
            "subject_smiles_dictionary_" + args.naming + ".yml",
        ),
        "w",
    ) as outfile:
        yaml.dump(smiles_dict_filelist, outfile, default_flow_style=False)

    # Write SMILES list to file
    clean_smiles_set = set(df_cleaned_smiles_patents.smiles)
    smiles_txt_file = os.path.join(args.output_dir, args.naming, f"smiles_list_{args.naming}.txt")
    with open(smiles_txt_file, "w") as f:
        for i, smi in enumerate(clean_smiles_set):
            if i == len(clean_smiles_set) - 1:
                f.write(smi)
            else:
                f.write(smi + "\n")

    return


if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args)
