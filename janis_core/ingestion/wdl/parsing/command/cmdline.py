
from dataclasses import dataclass, field
import regex as re 
import WDL

@dataclass
class CmdLine:
    elems: list[str | WDL.Expr.Placeholder] = field(default_factory=list)

    @property
    def is_set_command(self) -> bool:
        if isinstance(self.elems[0], str):
            if self.elems[0].startswith('set '):
                return True
        return False
    
    @property
    def is_shell_var_cmd(self) -> bool:
        if not isinstance(self.elems[0], str):
            return False
        if not re.match(r'^[ \t]*?(set|declare|env|export)', self.elems[0]):
            return False
        return True
    
    @property
    def is_comment(self) -> bool:
        if not isinstance(self.elems[0], str):
            return False
        if not re.match(r'^[ \t]*?#', self.elems[0]):
            return False
        return True