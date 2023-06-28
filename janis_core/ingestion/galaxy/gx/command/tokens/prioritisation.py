

from abc import ABC, abstractmethod
from .token import Token, TokenType


### ORDERING ###

class TokenOrderingStrategy(ABC):
    @abstractmethod
    def order(self, token_list: list[Token]) -> list[Token]:
        """selects a token in the list according to some prioritisation"""
        ...

class FirstTokenOrderingStrategy(TokenOrderingStrategy):
    def order(self, token_list: list[Token]) -> list[Token]:
        """selects a token in the list according to some prioritisation"""
        token_list.sort(key=lambda x: x.start)
        return token_list

class LongestTokenOrderingStrategy(TokenOrderingStrategy):
    def order(self, token_list: list[Token]) -> list[Token]:
        """selects a token in the list according to some prioritisation"""
        token_list.sort(key=lambda x: x.end - x.start)
        return token_list

class PriorityTokenOrderingStrategy(TokenOrderingStrategy):
    def order(self, token_list: list[Token]) -> list[Token]:
        priorities: dict[TokenType, int] = {
            TokenType.FUNCTION_CALL: 0,
            TokenType.BACKTICK_SHELL_STATEMENT: 0,
            TokenType.SCRIPT: 1,
            TokenType.GX_KW_DYNAMIC: 1,
            TokenType.GX_KW_STATIC: 2,
            TokenType.KV_LINKER: 3,
            TokenType.GX_INPUT: 4,
            TokenType.GX_OUTPUT: 4,
            TokenType.ENV_VAR: 5,
            TokenType.STRING: 8,
            TokenType.INTEGER: 9,
            TokenType.FLOAT: 9,
            TokenType.LINUX_TEE: 10,
            TokenType.LINUX_REDIRECT: 10,
            TokenType.LINUX_STREAM_MERGE: 10,
            TokenType.EMPTY_STRING: 11,
            TokenType.UNKNOWN: 999,
        }
        token_list.sort(key=lambda x: priorities[x.ttype])
        return token_list