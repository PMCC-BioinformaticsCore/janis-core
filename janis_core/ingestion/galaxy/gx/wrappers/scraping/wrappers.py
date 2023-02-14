

from abc import ABC, abstractmethod
import requests
import tarfile
import shutil
import os
from dataclasses import dataclass
from typing import Any, Optional
from xml.etree import ElementTree as et

from galaxy.tool_util.parser import get_tool_source
from janis_core.ingestion.galaxy.gx.wrappers import Wrapper
import gx.wrappers.scraping.utils as utils
from wrappers.WrapperCache import WrapperCache


THREADS = 10


def scrape_wrappers() -> None:
    scraper = WrapperScraper()
    scraper.scrape()


class ScrapeStrategy(ABC):
    def __init__(self, owner: str, repo: str, revision: str):
        self.owner = owner
        self.repo = repo
        self.revision = revision

    @abstractmethod
    def scrape(self) -> None:
        raise NotImplementedError()


@dataclass
class Revision:
    owner: str
    repo: str
    revision: str
    date_created: str
    download_folder: str = 'temp'

    @property
    def url(self) -> str:
        return f'https://toolshed.g2.bx.psu.edu/repos/{self.owner}/{self.repo}/archive/{self.revision}.tar.gz'

    def get_tar(self) -> tarfile.TarFile:
        response = requests.get(self.url, stream=True)
        tar = tarfile.open(fileobj=response.raw, mode='r:gz')
        tar.extractall(path=self.download_folder)
        return tar

    def get_wrappers(self) -> list[Wrapper]:
        tar = self.get_tar()
        folder = os.path.commonprefix(tar.getnames())
        xmlpaths = [f'{self.download_folder}/{x}' for x in tar.getnames() if x.endswith('.xml')]
        wrappers = self.extract_wrapper_info(xmlpaths)
        shutil.rmtree(f'{self.download_folder}/{folder}')
        return wrappers

    def extract_wrapper_info(self, xmlpaths: list[str]) -> list[Wrapper]:
        out: list[Wrapper] = []
        for path in xmlpaths:
            info = self.get_info(path)
            if info:
                info['owner'] = self.owner
                info['repo'] = self.repo
                info['revision'] = self.revision
                info['date_created'] = self.date_created
                info['requirements'] = []
                raise NotImplementedError()
                out.append(Wrapper(info))
        return out
    
    def get_info(self, path: str) -> Optional[dict[str, Any]]:
        tree = et.parse(path)
        root = tree.getroot()
        if root.tag == 'tool':
            tool_source = get_tool_source(path)  # TODO macro paths????
            return {
                'tool_id': tool_source.parse_id(),  # type: ignore
                'tool_build': tool_source.parse_version()  # type: ignore
            }
        return None



class WrapperScraper:

    def __init__(self):
        self.cache = WrapperCache()
        self.total_count: int = 0
        self.scraped_count: int = 0

    @property
    def percentage_complete(self) -> float:
        return (self.scraped_count / self.total_count) * 100

    def scrape(self) -> None:
        revisions = self.get_scrapable_revisions()
        self.total_count = len(revisions)
        print(f'revisions to scrape: {self.total_count}')
        for revision in revisions:
            try:
                self.scraped_count += 1
                print(f'\r{self.percentage_complete:0.1f}%', end='')
                wrappers = revision.get_wrappers()
                for wrapper in wrappers:
                    self.cache.add(wrapper)
            except Exception as e:
                print(e)

    def get_scrapable_revisions(self) -> list[Revision]:
        revision_data = utils.load_data(utils.REVISION_DATA_PATH)
        out: list[Revision] = []
        for owner, repo_data in revision_data.items():
            for repo, revision_list in repo_data.items():
                for revision in revision_list:
                    if not self.cache.get(revision=revision['revision']):
                        new_revision = Revision(
                            owner=owner,
                            repo=repo,
                            revision=revision['revision'],
                            date_created=revision['date_created'],
                        )
                        out.append(new_revision)
        return out
