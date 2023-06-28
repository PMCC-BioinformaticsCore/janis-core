
import os
import json
import xml.etree.ElementTree as et
from typing import Optional


def is_galaxy_workflow(path: str) -> bool:
    """checks that the provided path is a galaxy workflow"""
    with open(path, 'r') as fp:
        data = json.load(fp)
    if data['a_galaxy_workflow'] == 'true':
        return True
    return False

def is_tool_xml(path: str) -> bool:
    """checks that the provided path is tool xml"""
    root = et.parse(path).getroot()
    if root.tag == 'tool':
        return True
    return False

def is_macro_xml(path: str) -> bool:
    """checks that the provided path is tool xml"""
    root = et.parse(path).getroot()
    if root.tag == 'macros':
        return True
    return False

def get_xml_id(filepath: str) -> str:
    tree = et.parse(filepath)
    root = tree.getroot()
    return str(root.attrib['id']) # type: ignore

def get_xml_by_id(xmldir: str, query_id: str) -> Optional[str]:
    xmls = [x for x in os.listdir(xmldir) if x.endswith('.xml') and 'macros' not in x]
    for xml in xmls:
        filepath = f'{xmldir}/{xml}'
        tool_id = get_xml_id(filepath)
        if query_id == tool_id:
            return xml
    return None 

def get_macros(xmldir: str) -> list[str]:
    out: list[str] = []    
    xmls = [x for x in os.listdir(xmldir) if x.endswith('.xml')]
    for xml in xmls:
        path = f'{xmldir}/{xml}'
        if is_macro_xml(path):
            out.append(path)
    return out
