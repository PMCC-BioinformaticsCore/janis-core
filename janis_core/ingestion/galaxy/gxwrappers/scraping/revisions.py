

from datetime import datetime
from typing import Any, Tuple

from bs4 import Tag

from janis_core.ingestion.galaxy.runtime.dates import TOOLSHED_DATE_FMT, JANIS_DATE_FMT
import gx.wrappers.scraping.utils as utils

NEW_REPOS_ONLY = True
THREADS = 10

def scrape_revisions() -> None:
    scraper = RevisionScraper()
    scraper.scrape()

def get_revisions_in_repo(owner: str, repo: str) -> list[dict[str, str]]:
    URL = format_repo_url(owner, repo)
    page = utils.get_soup(URL)
    table = utils.get_table(page)
    ages = utils.get_table_col(table, 0)
    links = utils.get_table_col(table, 2)
    return extract_revisions_info(ages, links)

def format_repo_url(owner: str, repo: str) -> str:
    return f'https://toolshed.g2.bx.psu.edu/repos/{owner}/{repo}/'

def extract_revisions_info(ages: list[Tag], links: list[Tag]) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for age, link in zip(ages, links):
        revision: str = link.find('a').attrs['href'].rsplit('/', 1)[-1] # type: ignore
        date_created = datetime.strptime(age.text, TOOLSHED_DATE_FMT) # type: ignore
        date_string: str = date_created.strftime(JANIS_DATE_FMT)  # TODO does this mess up the time since no UTC offset?
        out.append({
            'revision': revision,
            'date_created': date_string,
        })
    return out


class RevisionScraper:
    """
    coordinates the scraping process
    handles file loading, updating saving
    """
    def __init__(self):
        self.local_data: dict[str, Any] = utils.load_data(utils.REVISION_DATA_PATH)
        self.total_count: int = 0
        self.scraped_count: int = 0

    @property
    def percentage_complete(self) -> float:
        return (self.scraped_count / self.total_count) * 100

    def scrape(self) -> None:
        owners_repos = self.get_scrapable_repos()
        self.total_count = len(owners_repos)
        print(f'repos to scrape: {self.total_count}')
        for owner, repo in owners_repos:
            self.try_scrape_repo_revisions(owner, repo)
            self.scraped_count += 1
            print(f'\r{self.percentage_complete:0.1f}%', end='')

    def get_scrapable_repos(self) -> list[Tuple[str, str]]:
        repo_data = utils.load_data(utils.REPO_DATA_PATH)
        out: list[Tuple[str, str]] = []
        for owner, repo_list in repo_data.items():
            for repo in repo_list:
                if self.should_scrape_repo_revisions(owner, repo):
                    out.append((owner, repo))
        return out

    def should_scrape_repo_revisions(self, owner: str, repo: str) -> bool:
        # if we only want new repos, don't scrape existing repos
        if NEW_REPOS_ONLY:
            if owner in self.local_data:
                if repo in self.local_data[owner]:
                    return False
        return True

    def try_scrape_repo_revisions(self, owner: str, repo: str) -> None:
        if owner not in self.local_data:
            self.local_data[owner] = {}
        try:
            self.local_data[owner][repo] = get_revisions_in_repo(owner, repo)  # takes time - do threaded
            utils.save_data(utils.REVISION_DATA_PATH, self.local_data) # save
        except Exception as e:
            print(e)

