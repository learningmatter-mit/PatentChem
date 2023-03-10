import argparse
import os
import tarfile
import urllib.request
from urllib.request import urlopen as uReq
from zipfile import ZipFile

from bs4 import BeautifulSoup as soup
from tqdm import tqdm


def create_directory(dir_path):
    if os.path.exists(dir_path):
        print(f"Directory for {dir_path} already exists.")
    else:
        os.makedirs(dir_path)
        print(f"Directory {dir_path} created")


def scrape_uspto(patent_year):
    # Define the USPTO web page according to the year:
    patent_url = os.path.join("https://bulkdata.uspto.gov/data/patent/grant/redbook/", patent_year)

    # Open a connection and download the webpage:
    uClient = uReq(patent_url)

    # Read the html of the page:
    page_html = uClient.read()

    # Close the webpage:
    uClient.close()

    # Parse the html using Beautiful Soup:
    page_soup = soup(page_html, "html.parser")

    # Get all the attributes of the page, containing the link to the weekly patent files
    patent_weekly_releases = page_soup.findAll("a")

    # Here, the link are selected from containers (ex: 'I20180102.tar').
    # They can be used to complete the url for download, e.g.
    # https://bulkdata.uspto.gov/data/patent/grant/redbook/2018/I20180102.tar
    patent_link = []
    for release in patent_weekly_releases:
        match = [s for s in release if ".tar" in s]
        matchzip = [s for s in release if ((".zip" in s) or (".ZIP" in s))]
        if match != []:
            if match[0].endswith(".tar"):
                patent_link.append(match)
        if matchzip != []:
            if (
                matchzip[0].lower().endswith(".zip")
                and ("SUPP" not in matchzip[0])
                and ("Grant" not in matchzip[0])
                and ("Red Book Viewer Public" not in matchzip[0])
            ):
                patent_link.append(matchzip)

    print(f"Found {len(patent_link)} releases from {patent_year}")

    return patent_link


class DownloadProgressBar(tqdm):
    # Progression bar for the download
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


def download_url(url, output_path):
    with DownloadProgressBar(unit="B", unit_scale=True, miniters=1, desc=url.split("/")[-1]) as t:
        urllib.request.urlretrieve(url, filename=output_path, reporthook=t.update_to)


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
        "--force_redownload",
        action="store_true",
        help="Use if you want to redownload files that already exist",
    )
    parser.add_argument(
        "--remove_compressed",
        action="store_true",
        help="Use if you want to remove the original .tar and .zip files for the weekly releases",
    )
    return parser


def main(args):
    # Change to directory that will contain all years
    os.chdir(args.data_dir)

    # Construct years list from parsed args
    if args.years == ["all"]:
        print("Preparing to download all USPTO patents from 2001 to 2023...")
        years = list(map(str, range(2001, 2024)))
    else:
        print("Preparing to download all USPTO patents from", ", ".join(args.years), "...")
        years = args.years

    # Iterate through years
    for patent_year in years:
        # Get list of weekly releases from web scraper
        patent_link = scrape_uspto(patent_year)

        # Create the directory for the year if it doesn't exist yet
        create_directory(patent_year)

        # Iterate through weekly releases
        for link in patent_link:
            year_link = os.path.join(patent_year, link[0])
            target_dir = year_link[:-4]

            # Create target directory to untar/unzip the files
            create_directory(target_dir)

            # Start downloading
            if not os.path.isfile(year_link) or args.force_redownload:
                download_url(
                    f"https://bulkdata.uspto.gov/data/patent/grant/redbook/{year_link}",
                    year_link,
                )
            else:
                print(
                    f"File {year_link} already exists. Use option '--force_redownload' if you want to download anyway."
                )

            # Once downloaded, untar/unzip in the corresponding directory and remove the original tar/zip file
            if link[0].endswith(".tar"):
                tar = tarfile.open(name=year_link, mode="r")
                tar.extractall(target_dir)
                tar.close()
            elif link[0].lower().endswith(".zip"):
                with ZipFile(year_link, "r") as zf:
                    zf.extractall(path=target_dir)

            if args.remove_compressed:
                os.remove(year_link)

    return


if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args)
