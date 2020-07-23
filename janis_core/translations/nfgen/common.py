from enum import Enum
from abc import ABC, abstractmethod
from typing import Optional


def filter_null(iterable):
    if iterable is None:
        return None
    if not isinstance(iterable, list):
        raise Exception(
            f"filter_null: Expected iterable ('{str(iterable)}') to be of type list, received '{type(iterable)}'"
        )

    return [el for el in iterable if el is not None]


def get_nf_string(value):
    if isinstance(value, list):
        return ", ".join(get_nf_string(v for v in value))

    elif isinstance(value, Enum):
        return value.value

    elif isinstance(value, NFBase):
        return value.get_string()

    return str(value)


class NFBase(ABC):
    @abstractmethod
    def get_string(self):
        raise Exception("Subclass must override .get_string() method")


class Channel(NFBase):
    def __init__(self, method, value):
        self.method = method
        self.value = value

    def get_string(self):
        value = get_nf_string(self.value)
        return f"Channel.{self.method}( {value} )"

    @classmethod
    def from_(cls, value):
        return cls("from", value)
