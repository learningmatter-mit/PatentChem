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
        description="Reset each given year's directory to its post-download, pre-select_chem state"
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
    parser.add_argument(
        "--remove_uncompressed_only",
        action="store_true",
        help="Use if you want to remove all uncompressed files but not uncompress the zip/tar files again",
    )
    return parser


def main(args):
    os.chdir(args.data_dir)
    for patent_year in args.years:
        os.chdir(patent_year)

        # make sure there are ZIP/zip/tar files present before deleting anything
        files = glob.glob("*")
        if not any([file.endswith(".tar") or file.lower().endswith(".zip") for file in files]):
            print(
                f"No zip/tar files found in {patent_year} directory. Reset is only possible if `--remove_compressed` was not used in download.py. Skipping year {patent_year}..."
            )
            continue

        # remove all directories but leave all ZIP/zip/tar files
        for file in glob.glob("*"):
            if os.path.isdir(file):
                os.system(f"rm -r {file}")

        # unzip/untar files for those years (assumes `--remove_compressed` was not used in download.py)
        if not args.remove_uncompressed_only:
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
