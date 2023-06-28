

import os
import json
import requests
from typing import Any
from bs4 import BeautifulSoup, Tag

# FILE IO 


def load_data(path: str) -> dict[str, Any]:
    """lockfile?"""
    if not os.path.exists(path):
        return {}
    with open(path, 'r') as fp:
        return json.load(fp)

def save_data(path: str, data: dict[str, Any]) -> None:
    """lockfile?"""
    with open(path, 'w') as fp:
        json.dump(data, fp)


# REQUESTS / HTML PARSING
def get_soup(url: str) -> BeautifulSoup:
    response = requests.get(url)
    return BeautifulSoup(response.content, "html.parser")

def get_table(page: BeautifulSoup) -> Tag:
    table = page.find(name='table')
    if not isinstance(table, Tag):
        raise RuntimeError()
    return table

def get_table_col(table: Tag, col: int) -> list[Tag]:
    rows = table.find_all('tr')
    rows = [row for row in rows if row.parent.name == 'tbody']
    elems = [row.find_all('td')[col] for row in rows]
    return elems

def format_iframe_url(iframe_elem: Tag) -> str:
    return f'https://toolshed.g2.bx.psu.edu/{iframe_elem.attrs["src"]}'


