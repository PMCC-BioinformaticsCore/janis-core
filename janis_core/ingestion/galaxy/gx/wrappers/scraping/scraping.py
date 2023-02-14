

"""

module info

main calls and operations here using ToolshedScraper etc

"""

import sys
sys.path.append('./src')

from janis_core.ingestion.galaxy.gx.wrappers.scraping.repositories import scrape_repos
from janis_core.ingestion.galaxy.gx.wrappers.scraping.revisions import scrape_revisions
from janis_core.ingestion.galaxy.gx.wrappers.scraping.wrappers import scrape_wrappers


def main(argv: list[str]):
    mode = argv[0]
    match mode:
        case 'repos':
            scrape_repos()
        case 'revisions':
            scrape_revisions()
        case 'wrappers':
            scrape_wrappers()
        case _:
            print(f'invalid mode {mode}')
            sys.exit()


if __name__ == '__main__':
    main(sys.argv[1:])

