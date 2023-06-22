

from bs4 import Tag
from collections import defaultdict
import gx.wrappers.scraping.utils as utils



def scrape_repos() -> None:
    """
    only single thread possible
    returns list of urls, where each url is a repo
    consists of all repos in the toolshed
    """
    URL = 'https://toolshed.g2.bx.psu.edu/repos'
    page = utils.get_soup(URL)
    table = utils.get_table(page)
    elems = utils.get_table_col(table, 0)
    repo_data = extract_repo_information(elems)
    utils.save_data(utils.REPO_DATA_PATH, repo_data)

def extract_repo_information(elems: list[Tag]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = defaultdict(list)
    for elem in elems:
        path = elem.find('a').attrs['href'].split('/') # type: ignore
        owner, repo = path[2], path[3]
        out[owner].append(repo) # type: ignore
    return out


