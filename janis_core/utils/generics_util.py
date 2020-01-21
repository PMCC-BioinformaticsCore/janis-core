import typing


is_py38 = hasattr(typing, "get_args")
is_py37 = not (is_py38) and hasattr(typing, "_GenericAlias")
is_py36 = not (is_py38 or is_py37) and hasattr(typing, "_Union")
is_py35 = not (is_py38 or is_py37 or is_py36)

if is_py38:

    def _is_generic(cls):
        if isinstance(cls, typing._GenericAlias):
            return True

        if isinstance(cls, typing._SpecialForm):
            return cls not in {typing.Any}

        return False

    def _is_base_generic(cls):
        args = typing.get_args(cls)
        args_are_empty = len(list(a for a in args if a)) == 0
        return args_are_empty or str(args[0]) == "~T"


elif is_py37:
    # python 3.7
    def _is_generic(cls):
        if isinstance(cls, typing._GenericAlias):
            return True

        if isinstance(cls, typing._SpecialForm):
            return cls not in {typing.Any}

        return False

    def _is_base_generic(cls):
        if isinstance(cls, typing._GenericAlias):
            if cls.__origin__ in {typing.Generic, typing._Protocol}:
                return False

            if isinstance(cls, typing._VariadicGenericAlias):
                return True

            return len(cls.__parameters__) > 0

        if isinstance(cls, typing._SpecialForm):
            return cls._name in {"ClassVar", "Union", "Optional"}

        return False


elif is_py36:
    # python 3.6
    def _is_generic(cls):
        return isinstance(
            cls, (typing.GenericMeta, typing._Union, typing._Optional, typing._ClassVar)
        )

    def _is_base_generic(cls):
        if isinstance(cls, (typing.GenericMeta, typing._Union)):
            return cls.__args__ in {None, ()}

        return isinstance(cls, typing._Optional)


else:
    # python 3.5
    def _is_generic(cls):
        return isinstance(
            cls,
            (
                typing.GenericMeta,
                typing.UnionMeta,
                typing.OptionalMeta,
                typing.CallableMeta,
                typing.TupleMeta,
            ),
        )

    def _is_base_generic(cls):
        if isinstance(cls, typing.GenericMeta):
            return all(isinstance(arg, typing.TypeVar) for arg in cls.__parameters__)

        if isinstance(cls, typing.UnionMeta):
            return cls.__union_params__ is None

        if isinstance(cls, typing.TupleMeta):
            return cls.__tuple_params__ is None

        if isinstance(cls, typing.CallableMeta):
            return cls.__args__ is None

        if isinstance(cls, typing.OptionalMeta):
            return True

        return False


def is_generic(cls):
    """
    Detects any kind of generic, for example `List` or `List[int]`. This includes "special" types like
    Union and Tuple - anything that's subscriptable, basically.
    """
    return _is_generic(cls)


def is_base_generic(cls):
    """
    Detects generic base classes, for example `List` (but not `List[int]`)
    """
    return _is_base_generic(cls)


def is_qualified_generic(cls):
    """
    Detects generics with arguments, for example `List[int]` (but not `List`)
    """
    return is_generic(cls) and not is_base_generic(cls)
