

from janis_core import CommandToolBuilder, WorkflowBuilder, InputNodeSelector, StepOutputSelector


def wrap_tool_in_workflow(tool: CommandToolBuilder) -> WorkflowBuilder:
    wf = WorkflowBuilder(f'{tool.id()}_wf')

    # add workflow input for each tool input
    for inp in tool._inputs:
        wf.input(
            identifier=inp.id(),
            datatype=inp.input_type,
            default=None,
        )

    # add step
    kwargs = {inp.id(): InputNodeSelector(wf.input_nodes[inp.id()]) for inp in tool._inputs}
    step = wf.step(
        f'{tool.id()}_step',
        tool(**kwargs)
    )

    # add output for each tool output
    for out in tool._outputs: 
        wf.output(
            identifier=out.id(),
            datatype=out.output_type,
            source=StepOutputSelector(step, out.id())
        )

    return wf 