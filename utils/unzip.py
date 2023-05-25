import argparse
import glob
import os
import tarfile
from zipfile import ZipFile


def create_directory(dir_path):
    if os.path.exists(dir_path):
        print(f"Directory for {dir_path} already exists.")
    else:
        os.makedirs(dir_path)
        print(f"Directory {dir_path} created")


def get_parser():
    parser = argparse.ArgumentParser(
        description="Unzips all zip/tar files for a given set of years. Can use after download.py if `--no_uncompress` was used."
    )
    parser.add_argument(
        "--years",
        type=str,
        nargs="+",
        required=True,
        help="Year(s) of patent files to unzip",
    )
    parser.add_argument(
        "--data_dir",
        type=str,
        default=".",
        help="Path where all data is be stored (e.g. /data/patents_data/)",
    )
    return parser


def main(args):
    os.chdir(args.data_dir)
    for patent_year in args.years:
        os.chdir(patent_year)

        # check if any non-zip/tar files are present
        files = glob.glob("*")
        if any([not (file.endswith(".tar") or file.lower().endswith(".zip")) for file in files]):
            print(
                f"Non-zip/tar files found in {patent_year} directory. Skipping year {patent_year}..."
            )
            continue

        for file in glob.glob("*"):
            target_dir = file[:-4]
            create_directory(target_dir)

            if file.endswith(".tar"):
                tar = tarfile.open(name=file, mode="r")
                tar.extractall(target_dir)
                tar.close()
            elif file.lower().endswith(".zip"):
                with ZipFile(file, "r") as zf:
                    zf.extractall(path=target_dir)
    return


if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args)
