


import json
import requests
from typing import Any, Optional
from requests import Response

from janis_core.ingestion.galaxy.logs import logging



def make_api_request(request_url: str) -> Any:
    response = make_request(request_url)
    return handle_response(response)

def make_request(request_url: str) -> Optional[Response]:
    # make requests to get information about tools with similar name
    iterations = 1
    response: Optional[Response] = None
    while response is None:
        if iterations > 5:
            break
        try:
            response = requests.get(request_url, timeout=5)
        except requests.exceptions.Timeout:
            response = None
        iterations += 1
    return response

def handle_response(response: Optional[Response]) -> Any:
    if response is None or response.status_code != 200:
        logging.no_ga4gh_data()
        return None
    else:
        return json.loads(response.text)