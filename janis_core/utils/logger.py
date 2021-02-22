"""
    Logger - Controls logging of the application
"""
import sys
import os
from datetime import datetime
from typing import Optional, TextIO, Dict, Tuple
import traceback


class _bcolors:
    """
    Lookup table: https://en.wikipedia.org/wiki/ANSI_escape_code#3/4_bit
    """

    HEADER = "\033[95m"  # Bright magenta
    OKBLUE = "\033[94m"  # Bright blue
    OKGREEN = "\033[92m"  # Bright green
    WARNING = "\033[93m"  # Bright yellow
    FAIL = "\033[91m"  # Bright red
    ENDC = "\033[0m"  # SGR (Reset / Normal)
    BOLD = "\033[1m"  # SGR (Bold or increased intensity
    ITALIC = "\033[3m"  # SGR (Italic)
    UNDERLINE = "\033[4m"  # SGR (Underline)


class LogLevel:
    CRITICAL = 2  # RED
    WARNING = 3  # YELLOW
    INFO = 4  # WHITE
    DEBUG = 5  # GREY
    VERBOSE = 6  # GREY

    @staticmethod
    def get_color(level: int):
        if level == LogLevel.CRITICAL:
            return _bcolors.FAIL

        if level == LogLevel.WARNING:
            return _bcolors.WARNING

        if level == LogLevel.INFO:
            return _bcolors.OKBLUE

        if level == LogLevel.DEBUG:
            return _bcolors.ENDC

        if level == LogLevel.VERBOSE:
            return _bcolors.ENDC

        return _bcolors.ENDC

    @staticmethod
    def get_str(level: int):
        if level == LogLevel.CRITICAL:
            return "CRITICAL"

        if level == LogLevel.WARNING:
            return "WARN"

        if level == LogLevel.INFO:
            return "INFO"

        if level == LogLevel.DEBUG:
            return "DEBUG"

        if level == LogLevel.VERBOSE:
            return "VERB"

        return ""

    @staticmethod
    def from_str(level: str):
        level = level.upper()
        if level == "CRITICAL":
            return LogLevel.CRITICAL
        if level == "WARN":
            return LogLevel.WARNING
        if level == "INFO":
            return LogLevel.INFO
        if level == "DEBUG":
            return LogLevel.DEBUG
        if level.startswith("VERB"):
            return LogLevel.VERBOSE
        if level == "NONE":
            return None

        raise KeyError(
            f"Couldn't detect LogLevel type '{level}', expected one of: "
            f"CRITICAL, WARN, INFO, DEBUG, NONE"
        )


class Logger:
    CONSOLE_LEVEL: Optional[int] = LogLevel.INFO
    __TEMP_CONSOLE_LEVEL: Optional[int] = None

    WRITE_LEVELS: Dict[int, Tuple[Optional[str], Optional[TextIO]]] = {}

    # DEPRECATED
    WRITE_LEVEL: Optional[int] = LogLevel.DEBUG

    WRITE_LOCATION: Optional[str] = None
    __WRITE_POINTER: Optional[TextIO] = None

    last_write: datetime = datetime.now()

    @staticmethod
    def set_console_level(level: Optional[int]):
        Logger.CONSOLE_LEVEL = level

    @staticmethod
    def mute():
        """
        Decreases log level to critical until unmute is called
        """
        Logger.__TEMP_CONSOLE_LEVEL = Logger.CONSOLE_LEVEL
        Logger.set_console_level(LogLevel.CRITICAL)

    @staticmethod
    def unmute():
        """
        unmute the console, if __TEMP is None, should not do anything
        :return:
        """
        if Logger.__TEMP_CONSOLE_LEVEL is not None:
            Logger.set_console_level(Logger.__TEMP_CONSOLE_LEVEL)

    # Deprecated
    @staticmethod
    def set_write_level(level: Optional[int]):
        Logger.WRITE_LEVELS = {level: (Logger.WRITE_LOCATION, Logger.__WRITE_POINTER)}
        Logger.WRITE_LEVEL = level

    @staticmethod
    def set_write_location(location: str):
        for p in Logger.WRITE_LEVELS.values():
            (path, pointer) = p
            if pointer is not None and not pointer.closed:
                pointer.close()

        if location is None:
            Logger.WRITE_LEVELS = {}
        else:
            level = Logger.WRITE_LEVEL
            # legacy
            Logger.WRITE_LOCATION = location
            Logger.__WRITE_POINTER = open(location, "a")
            Logger.WRITE_LEVELS = {level: (location, Logger.__WRITE_POINTER)}

    @staticmethod
    def close_file():
        nwl = {}
        for k, p in Logger.WRITE_LEVELS.items():
            path, pointer = p
            nwl[k] = (path, None)
            if pointer is not None and not pointer.closed:
                pointer.close()

        Logger.WRITE_LEVELS = nwl
        Logger.__WRITE_POINTER = None

    @staticmethod
    def log(message: str, level: int = LogLevel.VERBOSE):
        if level is None:
            # This is a developer error, we should never try to log with no level, it's purely for
            return

        m = f"{Logger.get_prefix(level)}: {message}"
        if Logger.CONSOLE_LEVEL is not None and level <= Logger.CONSOLE_LEVEL:
            print(LogLevel.get_color(level) + m + _bcolors.ENDC, file=sys.stderr)

        # if level <= LogLevel.CRITICAL:
        #     traceback.print_stack(limit=12)
        #     raise Exception(traceback.extract_stack(limit=5))

        should_write = (datetime.now() - Logger.last_write).total_seconds() >= 1

        for loglevel, p in Logger.WRITE_LEVELS.items():
            pointer = p[1]
            if level > loglevel or pointer is None or pointer.closed:
                continue
            pointer.write(m + "\n")

            if should_write:
                Logger.last_write = datetime.now()
                pointer.flush()
                os.fsync(pointer.fileno())

    @staticmethod
    def debug(message: str):
        Logger.log(message, LogLevel.DEBUG)

    @staticmethod
    def info(message: str):
        Logger.log(message, LogLevel.INFO)

    @staticmethod
    def warn(message: str):
        Logger.log(message, LogLevel.WARNING)

    @staticmethod
    def critical(message: str):
        Logger.log(message, LogLevel.CRITICAL)

    @staticmethod
    def log_ex(ex: Exception):
        Logger.critical(repr(ex))

    @staticmethod
    def guess_log(message: str, default_level=LogLevel.WARNING):
        if message.startswith("DEBUG"):
            Logger.debug(message[min(len(message) - 1, 6) :])
        elif message.startswith("INFO"):
            Logger.info(message[min(len(message) - 1, 5) :])
        elif message.startswith("WARNING"):
            Logger.warn(message[min(len(message) - 1), 8:])
        elif message.startswith("CRITICAL"):
            Logger.critical(message[min(len(message) - 9), 6:])

        else:
            Logger.log(message, level=default_level)

    @staticmethod
    def get_prefix(level: int):
        return f"{datetime.now().replace(microsecond=0).isoformat()} [{LogLevel.get_str(level)}]"
