from abc import ABC


class Selector(ABC):
    pass


class WildcardSelector(Selector):
    def __init__(self, wildcard, first_element=None):
        self.wildcard = wildcard
        self.first_element = first_element
