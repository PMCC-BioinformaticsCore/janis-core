

from textwrap import indent


NF_INDENT = '  '

def call_fmt0(name: str) -> str:
    return f'{name}()\n'

def call_fmt1(name: str, input: str) -> str:
    return f'{name}( {input} )\n'

def call_fmt2(name: str, inputs: list[str]) -> str:
    call_str = f'{name}(\n'
    for i, inp in enumerate(inputs):
        comma = ',' if i < len(inputs) - 1 else ''
        call_str += f'{NF_INDENT}{inp}{comma}\n'
    call_str += ')\n'
    return call_str

def format_process_call(name: str, inputs: list[str], ind: int=0) -> str:
    if len(inputs) == 0:
        call_str = call_fmt0(name)
    elif len(inputs) == 1:
        call_str = call_fmt1(name, inputs[0])
    elif len(inputs) > 1:
        call_str = call_fmt2(name, inputs)
    else:
        raise RuntimeError()

    return indent(call_str, ind * NF_INDENT)





