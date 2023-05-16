




def looks_like_regex(pattern: str) -> bool:
    """Return True if the pattern looks like a regex"""
    return any(char in pattern for char in '*?[]{}()|')