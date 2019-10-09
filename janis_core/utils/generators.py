from uuid import uuid4

scattered_ordered_variable_identifiers = [
    "i",
    "j",
    "k",
    "x",
    "y",
    "z",
    "a",
    "b",
    "c",
    "ii",
    "jj",
    "kk",
    "xx",
    "yy",
    "zz",
]


def generate_new_id_from(source, forbidden_ids):
    new_var = source[0]
    max_size_of_subid = None
    while new_var in forbidden_ids:
        if len(scattered_ordered_variable_identifiers) > 0:
            new_var = scattered_ordered_variable_identifiers.pop(0)
        elif max_size_of_subid == len(source):
            # generate guid
            new_var = str(uuid4()).replace("-", "_")
        else:
            max_size_of_subid = (max_size_of_subid or 0) + 1
            new_var = source[-max_size_of_subid:]
    return new_var
