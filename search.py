import argparse
import os
import pickle
import xml.etree.ElementTree as ET

from rdkit import Chem


def directory_creation(dir_path):
    try:
        # Create target Directory
        os.makedirs(dir_path)
        print("Directory ", dir_path, " Created ")
    except FileExistsError:
        print("Directory ", dir_path, " already exists")


def file_open(patent_year, data_dir):

    # To open the file
    with open(os.path.join(data_dir, "chemistry_patent_list_" + patent_year), "rb") as fp:
        chemical_patent_list = pickle.load(fp)

    with open(os.path.join(data_dir, "smile_normal_list_" + patent_year), "rb") as fp:
        smile_normal_list = pickle.load(fp)

    with open(os.path.join(data_dir, "path_mol_normal_list_" + patent_year), "rb") as fp:
        path_mol_normal_list = pickle.load(fp)

    with open(os.path.join(data_dir, "path_mol_special_list_" + patent_year), "rb") as fp:
        path_mol_special_list = pickle.load(fp)

    with open(os.path.join(data_dir, "mol_special_list_" + patent_year), "rb") as fp:
        mol_special_list = pickle.load(fp)

    return (
        chemical_patent_list,
        smile_normal_list,
        path_mol_normal_list,
        path_mol_special_list,
        mol_special_list,
    )


def parsingXML(
    subject_list,
    data_dir,
    chemical_patent_list,
):
    print("Subject list:", subject_list)
    # Parsing the XML file containing the claim of the patent
    counter = 0
    index_list = []
    code_list = []
    for patent in chemical_patent_list:

        patent_ = os.path.join(data_dir, patent)

        try:
            try:
                claim = ET.ElementTree(file=patent_ + patent_[-20:-1] + ".XML")
            except FileNotFoundError as e:
                print(e)
                pass
        except ET.ParseError:
            print(f"Failed with ET.ParseError on {patent_}")
            continue

        # For looking into only the abstract
        abstract_xml = claim.find("abstract")
        abstract_text = ET.tostring(abstract_xml, method="text").decode("utf-8")

        match = 0
        for subject in subject_list:
            if subject in abstract_text:
                match += 1

        # Build a list (xml_text) with the text contained in the title and bulk text of the xml file.
        xml_text = []
        for tag_list in ["invention-title", "p"]:
            for elem in claim.iter(tag=tag_list):
                xml_text.append(elem.text)

        # Check for words matching the target (subject_list)
        for text in xml_text:
            if text is not None:
                for subject in subject_list:
                    if subject in text:
                        match += 1

        if match > 1:  # If match, add the index of the patent to the list (index_list)
            counter += 1
            # Get the patent code
            for elem in claim.iter(tag="country"):
                country = elem.text
                break  # The right one is on the first iteration
            for elem in claim.iter(tag="doc-number"):
                docnumber = elem.text
                break
            for elem in claim.iter(tag="kind"):
                kind = elem.text
                break
            if docnumber.startswith(
                "0"
            ):  # Again, consistency! docnumber = 09862000 => patent code = US9862000B2
                patent_code = country + docnumber[1:] + kind
            else:
                patent_code = country + docnumber + kind

            # Index list
            try:
                index_list.append(chemical_patent_list.index(patent))
            except ValueError as e:
                print("ValueError:", e)
                pass
            # Code list
            code_list.append(patent_code)

    return index_list, code_list


def dictionary(
    chemical_patent_list,
    smile_normal_list,
    path_mol_normal_list,
    path_mol_special_list,
    mol_special_list,
    index_list,
    code_list,
):
    # Extract all molecules (path) of the patent from the patent path
    index_mol_list = []
    code_mol_list = (
        []
    )  # When we extract the molecules, this is a cheap trick to keep the code patent
    for item in range(len(index_list)):
        extracted_mol = [
            i
            for i, x in enumerate(path_mol_normal_list)
            if chemical_patent_list[index_list[item]] in x
        ]
        # Done with intermediary var To keep the length
        index_mol_list.append(extracted_mol)
        for j in range(len(extracted_mol)):
            code_mol_list.append(code_list[item])

    index_mol_flat_list = [item for sublist in index_mol_list for item in sublist]

    # Building of the smile list from the correlated index of path_mol_normal_list and smile_normal_list
    subject_smiles_list = []
    for index in index_mol_flat_list:
        subject_smiles_list.append(smile_normal_list[index])

    """BUILDING THE DICTIONARY - key: the smiles, property: the patent code

    1) Currently, the patent code and the smiles are correlated by the index of their respective lists
    "subject_smiles_list" and "code_list"

    2) We need to split the Smiles that contains '.'. Within the loop, the smile has to be reprocessed
    with RDKit and be added to the dictionary with its corresponding code (that iterates with the for loop)

    """

    subject_smiles_dictionary = {}
    for index in range(len(subject_smiles_list)):
        # split clusters (SMILES separated by '.')
        subject_smiles_splitted = subject_smiles_list[index].split(".")
        for smile in subject_smiles_splitted:
            # Reprocess with RDKit
            resmile = Chem.MolFromSmiles(smile)
            if resmile is None:
                continue
            else:
                resmile_processed = Chem.MolToSmiles(resmile)

            subject_smiles_dictionary.setdefault(resmile_processed, []).append(code_mol_list[index])

            # if resmile_processed in subject_smiles_dictionary:
            #    subject_smiles_dictionary.setdefault(resmile_processed, []).append(code_mol_list[index])
            # else:
            #    subject_smiles_dictionary[resmile_processed] = code_mol_list[index]

    return subject_smiles_dictionary


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
        "--data_dir",
        type=str,
        default=".",
        help="Path where all data is be stored (e.g. /data/patents_data/)",
    )
    parser.add_argument(
        "--naming",
        type=str,
        required=True,
        help="Name of query (location where results will be stored)",
    )
    parser.add_argument(
        "--subject_list",
        type=str,
        nargs="+",
        required=True,
        help="Terms to query (include irregular plurals and variations on words)",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default=".",
        help="Parent directory of 'args.naming'",
    )
    return parser


def main(args):

    directory_creation(os.path.join(args.output_dir, args.naming))

    if args.years == ["all"]:
        print("Preparing to search chemistry patents from 2005 to 2023...")
        years = list(map(str, range(2005, 2023)))
        # TODO: doesn't currently work for 2001-2004 because of different XML structure
    else:
        print("Preparing to search chemistry patents from", ", ".join(args.years), "...")
        years = args.years

    for patent_year in years:
        if not os.path.exists(os.path.join(args.data_dir, patent_year)):
            print(f"Directory for year {patent_year} does not exist in this data_dir")
            continue

        if os.path.isfile(
            args.naming + "_".join(["/subject_smiles_dictionary", args.naming, patent_year])
        ):
            print("Year " + patent_year + " already processed.")
            continue
        else:
            print("Starting the search for " + patent_year)
            # File opening
            print("File opening")
            (
                chemical_patent_list,
                smile_normal_list,
                path_mol_normal_list,
                path_mol_special_list,
                mol_special_list,
            ) = file_open(patent_year, args.data_dir)

            # Parsing
            print("XML parsing")
            index_list, code_list = parsingXML(
                args.subject_list,
                args.data_dir,
                chemical_patent_list,
                smile_normal_list,
                path_mol_normal_list,
                path_mol_special_list,
                mol_special_list,
                patent_year,
            )
            print("index_list:", index_list)
            # Writing of the dictionary
            print("Writing of the dictionary")
            subject_smiles_dictionary = dictionary(
                chemical_patent_list,
                smile_normal_list,
                path_mol_normal_list,
                path_mol_special_list,
                mol_special_list,
                index_list,
                code_list,
            )

            with open(
                args.naming + "_".join(["/subject_smiles_dictionary", args.naming, patent_year]),
                "wb",
            ) as file_subject_smiles_dictionary:
                pickle.dump(subject_smiles_dictionary, file_subject_smiles_dictionary)
            print("len(subject_smiles_dictionary):", len(subject_smiles_dictionary))

    return


if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args)
