import re


class Validators:
    identifier_regex = r"[a-zA-Z][a-zA-Z0-9_]+\Z"
    compiled_identifier = re.compile(identifier_regex)

    extra_prohibited_keys = {
        "identifier",
        "tool",
        "scatter",
        "ignore_missing",
        "output",
        "input",
    }

    @staticmethod
    def validate_identifier(identifier) -> bool:
        if identifier in Validators.extra_prohibited_keys:
            return False
        a = Validators.compiled_identifier.match(identifier)
        return a is not None

    @staticmethod
    def reason_for_failure(identifier) -> str:
        if identifier in Validators.extra_prohibited_keys:
            return "it is a reserved keyword"
        a = Validators.compiled_identifier.match(identifier)
        if a is not None:
            return (
                f"it was not validated by '{Validators.identifier_regex}' "
                f"(must start with letters, and then only contain letters, numbers and an underscore)"
            )

        return "Undefined"
