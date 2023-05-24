import argparse
import os
import pickle
import shutil
from zipfile import ZipFile

from rdkit import Chem


def patent_directory(patent_year, patent_link, data_dir):
    """
    Define the path for the directories ("subdirectory")
    and the patents (contained in zip files, "list_zip")
    """
    print("Beginning Step 1: Path definition")

    list_path = []
    list_zip = []
    for link in patent_link:
        link_dir = os.path.join(data_dir, patent_year, link)
        link_link_ice_dir = os.path.join(link_dir, link, "project/pdds/ICEwithdraw/")
        link_ice_dir = os.path.join(link_dir, "project/pdds/ICEwithdraw/")

        # The following cases deal with the lack on consistency in the
        # file/directories structure of the bulk data files
        if os.path.isdir(
            os.path.join(link_dir, link, "DESIGN")
        ):  # Already in standardized structure
            current_path = os.path.join(link_dir, link)
        elif os.path.isdir(link_link_ice_dir):
            for dir in os.listdir(
                os.path.join(
                    link_link_ice_dir,
                    link,
                )
            ):  # remove unnecessary directories in 2008-2010
                shutil.move(
                    os.path.join(
                        link_link_ice_dir,
                        link,
                        dir,
                    ),
                    os.path.join(link_dir, link, dir),
                )
            shutil.rmtree(os.path.join(link_dir, link, "project"))
            current_path = link_dir
        elif os.path.isdir(link_ice_dir):
            for dir in os.listdir(
                os.path.join(
                    link_ice_dir,
                    link,
                )
            ):  # remove unnecessary directories in 2008-2010
                shutil.move(
                    os.path.join(
                        link_ice_dir,
                        link,
                        dir,
                    ),
                    os.path.join(link_dir, link, dir),
                )
            shutil.rmtree(os.path.join(link_dir, "project"))
            current_path = link_dir
        elif link == "I20201006_ST26_Sample":  # edge case
            current_path = os.path.join(data_dir, patent_year, link, "20201006")
        else:
            print(f"Unknown directory structure encountered for {link}")

        list_path.append(current_path)
        subdirectory = os.listdir(current_path)
        subdirectory = [
            s for s in subdirectory if ".txt" not in s
        ]  # don't consider any present .txt files
        # subdirectory_zip = [s for s in subdirectory if ".zip" in s]

        for item in subdirectory:
            subdirectory_zip = os.listdir(os.path.join(current_path, item))
            for element_zip in subdirectory_zip:
                if element_zip.lower().endswith(".zip"):
                    list_zip.append(os.path.join(current_path, item, element_zip))

    print("Step 1 Complete")

    return list_path, list_zip


# Unzipping the relevant files
def patent_unzipper(list_path, list_zip, remove_compressed, data_dir):
    # Recording of the patent names related with chemistry in a txt file
    chemistry_patent_list = []
    molfile_list = []
    mol_list = []
    # count = 0

    print("Beginning Step 2: patent selection")

    # 2.1) Loop through all zip files
    if len(list_zip) > 0:  # TODO: handle case where some zip files already parsed
        for i, currentzip in enumerate(
            list_zip
        ):  # currentzip contains the full path of every patent
            if not (i % 40000):
                frac = i / len(list_zip)
                print(f"Fraction of ZIP files complete: {frac:.2f}")  # status update
            if os.path.isfile(currentzip):
                # 2.1.1) Open zip and return the list of files contained (.XML/.CDX/.MOL/.TIF)
                with ZipFile(currentzip, "r") as zf:
                    zf_list = zf.namelist()

                    # 2.1.2) Select the patent that are related to chemistry (contains .CDX file)
                    test_cdx = False  # Reinitialize the tester
                    for k in range(1, len(zf_list)):
                        if zf_list[k].lower().endswith(".cdx"):
                            test_cdx = True
                            break

                    if test_cdx:  # extract the zip if chemistry related
                        # Register this patent in the list
                        current_name = currentzip[:-4] + currentzip[-24:-4] + "/"
                        current_name = current_name.replace(data_dir, "")
                        chemistry_patent_list.append(current_name)
                        if not os.path.isdir(currentzip[:-4] + "/"):
                            zf.extractall(path=currentzip[:-4] + "/")

                # Finish the loop with currentzip by deleteting the zip file.
                # The ones that didn't contain any .CDX file were not extracted.
                if remove_compressed:
                    os.remove(currentzip)

    print("Step 2 Complete")
    # 2.2) Check the extracted files (their path is recorded in chemistry_patent_list),
    # keep the .MOL/.XML/.CDX and remove .TIF.
    print("Beginning Step 3: extracted file parsing and cleaning")
    for patent in chemistry_patent_list:
        patent_ = os.path.join(data_dir, patent)
        for element in os.listdir(patent_):
            # smile_sublist = []
            if element.lower().endswith(
                ".mol"
            ):  # Extraction of the Mol using RDkit and saving in list + .txt
                # Record the path
                current_name = patent_ + element
                current_name = current_name.replace(data_dir, "")
                molfile_list.append(current_name)
                # Intermediate state, extract the molecule from the .MOL
                mol_list.append(Chem.MolFromMolFile(patent_ + element))

            # Remove the TIF files to save memory.
            elif element.lower().endswith(".tif"):  # Remove .TIF
                os.remove(patent_ + element)  # Comment this elif to preserve .TIF files
    print("Step 3 Complete")
    return chemistry_patent_list, molfile_list, mol_list


"""
- Extracting the SMILE of the molecules once and for all -

The SMILES will be listed all together and indexed.
When using the search engine, the SMILES will be collected from the SMILES python list
that result from these variables
(the cells below called SMILE writing are saving it into files with pickle)
"""


def patent_smile(chemistry_patent_list, molfile_list, mol_list, data_dir):
    # 3) Differentiation between normal molecules and special cases:
    # Special case => Failure of the extraction, may contain R or other special element.
    # Two layers of specificity were determined: with / without RDKit sanitization in Chem.MolFromMolFile
    print("Beginning Step 4: SMILES writing")

    molfile_special_list = []
    molfile_normal_list = []
    # smile_list = []
    smile_normal_list = []
    mol_normal_list = []
    mol_special_list = []
    bugcounter = 0

    for index in range(len(mol_list)):
        if mol_list[index] is None:
            mol_special_list.append(Chem.MolFromMolFile(mol_list[index], sanitize=False))
            current_name = molfile_list[index]
            current_name = current_name.replace(data_dir, "")
            molfile_special_list.append(current_name)
        else:
            mol_normal_list.append(mol_list[index])
            current_name = molfile_list[index]
            current_name = current_name.replace(data_dir, "")
            molfile_normal_list.append(current_name)

            # 2.4) Transform .MOL into SMILE and write the corresponding .txt file /
            # The index is related to the chemical_mol_list, which gives the path
            try:
                smile = Chem.MolToSmiles(mol_list[index])
            except Exception:
                bugcounter += 1
                smile = " "

            smile_normal_list.append(smile)
            with open(molfile_list[index][:-4] + ".txt", "w") as smile_file:
                smile_file.write(smile)

    # Location of fail
    # mol_fail_location = [i for i,x in enumerate(mol_list) if x == None]

    print("Step 4 Complete")
    return (
        molfile_special_list,
        molfile_normal_list,
        mol_normal_list,
        mol_special_list,
        smile_normal_list,
        bugcounter,
    )


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
        "--remove_compressed",
        action="store_true",
        help="Use if you want to remove all original .zip files after unzipping those relevant to chemistry",
    )
    return parser


def main(args):

    if args.years == ["all"]:
        print("Preparing to select chemistry patents from 2001 to 2023...")
        years = list(map(str, range(2001, 2024)))
    else:
        print("Preparing to select chemistry patents from", ", ".join(args.years), "...")
        years = args.years

    for patent_year in years:
        print(f"Starting year {patent_year}")
        patent_link = os.listdir(os.path.join(args.data_dir, patent_year))

        # Remove unwanted files from list
        if ".ipynb_checkpoints" in patent_link:
            patent_link.remove(".ipynb_checkpoints")
        if "nohup.out" in patent_link:
            patent_link.remove("nohup.out")
        patent_link = [
            f for f in patent_link if not (f.lower().endswith(".zip") or f.lower().endswith(".tar"))
        ]

        list_path, list_zip = patent_directory(patent_year, patent_link, args.data_dir)

        print("len(list_path):", len(list_path))
        print("len(list_zip):", len(list_zip))

        # Unzipping
        chemistry_patent_list, molfile_list, mol_list = patent_unzipper(
            list_path, list_zip, args.remove_compressed, args.data_dir
        )

        # List saving into binary files
        os.chdir(args.data_dir)

        # To record the lists into a file
        with open("chemistry_patent_list_" + patent_year, "wb") as file_chemistry_patent_list:
            pickle.dump(chemistry_patent_list, file_chemistry_patent_list)
        print("len(chemistry_patent_list):", len(chemistry_patent_list))

        with open("molfile_list_" + patent_year, "wb") as file_mol_patent_list:
            pickle.dump(molfile_list, file_mol_patent_list)
        print("len(molfile_list):", len(molfile_list))

        # SMILE writing
        print("len(chemistry_patent_list):", len(chemistry_patent_list))
        print("len(mol_list):", len(mol_list))

        (
            molfile_special_list,
            molfile_normal_list,
            mol_normal_list,
            mol_special_list,
            smile_normal_list,
            bugcounter,
        ) = patent_smile(chemistry_patent_list, molfile_list, mol_list, args.data_dir)

        # To record the lists into a file
        with open("path_mol_special_list_" + patent_year, "wb") as file_molfile_special_list:
            pickle.dump(molfile_special_list, file_molfile_special_list)
        print("len(molfile_special_list):", len(molfile_special_list))

        with open("path_mol_normal_list_" + patent_year, "wb") as file_molfile_normal_list:
            pickle.dump(molfile_normal_list, file_molfile_normal_list)
        print("len(molfile_normal_list):", len(molfile_normal_list))

        with open("mol_normal_list_" + patent_year, "wb") as file_mol_normal_list:
            pickle.dump(mol_normal_list, file_mol_normal_list)
        print("len(mol_normal_list):", len(mol_normal_list))

        with open("smile_normal_list_" + patent_year, "wb") as file_smile_normal_list:
            pickle.dump(smile_normal_list, file_smile_normal_list)
        print("len(smile_normal_list):", len(smile_normal_list))

        with open("mol_special_list_" + patent_year, "wb") as file_mol_special_list:
            pickle.dump(mol_special_list, file_mol_special_list)
        print("len(mol_special_list):", len(mol_special_list))

        print(f"Completed year {patent_year}")

    return


if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args)
