




class ProcessOutput:
    def __init__(
        self,
        qualifier: Qualifier,
        name: str,
        expression: str,
        is_optional: bool = False,
        into: Optional[Union[str, List[str]]] = None,
        attributes: Optional[Union[str, List[str]]] = None,
    ):

        self.qualifier = qualifier
        self.name = name
        self.expression = expression
        self.into = into
        self.optional = is_optional
        self.attributes = attributes

    def get_string(self) -> str:
        if self.qualifier != OutputProcessQualifier.tuple:
            if not utils.is_simple_path(self.expression):
                self.expression = f'"${{{self.expression}}}"'

        els = [self.qualifier.value, f"{self.expression}"]

        if self.optional is True:
            els.extend(["optional", "true"])

        if self.into:
            intochannels = (
                ", ".join(str(i) for i in self.into)
                if isinstance(self.into, list)
                else str(self.into)
            )
            els.extend(["into", intochannels])
        if self.attributes:
            if isinstance(self.attributes, list):
                els.extend(str(a) for a in self.attributes)
            else:
                els.append(str(self.attributes))

        els.append(f", emit: {self.name}")

        return " ".join(str(e) for e in els).strip()


class TupleElementForOutput:
    def __init__(
        self,
        qualifier: OutputProcessQualifier,
        expression: str,
    ):

        self.qualifier = qualifier
        self.expression = expression

    def get_string(self) -> str:
        return f'{self.qualifier.value}("${{{self.expression}}}")'



