

workflow_path: str
dev_partial_eval: bool

def set_path(value: str) -> None:
    global workflow_path
    workflow_path = value

def set_dev_partial_eval(value: bool) -> None:
    global dev_partial_eval
    dev_partial_eval = value

