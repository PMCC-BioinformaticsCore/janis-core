

from typing import Any
from janis_core.tool.commandtool import (
    CommandTool, 
    ToolInput, 
    ToolArgument
)

from .unwrap import unwrap_expression
from janis_core.types import (
    Filename,
    File,
    Array,
    Boolean
)



def gen_prescript_for_cmdtool(
    tool: CommandTool, 
    inputs: list[ToolInput], 
    resource_vars: list[str], 
    input_in_selectors: dict[str, Any]
) -> str:
    return ProcessPreScriptGenerator(
        tool,
        inputs,
        resource_vars,
        input_in_selectors,
    ).generate()



# main class 

class ProcessPreScriptGenerator:
    def __init__(
        self,
        tool: CommandTool, 
        inputs: list[ToolInput], 
        resource_vars: list[str], 
        input_in_selectors: dict[str, Any]
    ):
        assert(tool)
        self.tool = tool
        self.inputs = inputs        # exposed_inputs vs internal_inputs???
        self.resource_vars = resource_vars
        self.input_in_selectors = input_in_selectors

    def generate(self) -> str:
        lines: list[str] = []
        lines += self.gen_inputs_in_selector()
        lines += self.gen_expression_inputs()
        lines += self.gen_input_vars()
        return '\n'.join(lines)

    def gen_inputs_in_selector(self) -> list[str]:
        """
        If there is any input being referenced by Janis InputSelector,
        we need to add their Groovy variable definitions
        Grace comment: I'm not sure how this function works or why its needed. 
        """
        if self.tool.versioned_id() not in self.input_in_selectors:
            return []

        lines: list[str] = []
        input_keys = [i.id() for i in self.inputs]
        for k in self.input_in_selectors[self.tool.versioned_id()]:
            if k not in input_keys and k not in self.resource_vars:
                val = "''"
                code = f'def {k} = {val}'
                lines.append(code)
        return lines

    def gen_expression_inputs(self) -> list[str]:
        """
        Generate Groovy code to represent the values of input variable definitions for complex expressions
        """
        lines: list[str] = []
        for i in self.inputs:
            if hasattr(i, 'input_type'):
                input_type = i.input_type
            elif hasattr(i, 'intype'):
                input_type = i.intype
            else:
                raise Exception('Failed to get input type attribute')

            if isinstance(input_type, Filename):
                val = unwrap_expression(
                    value=input_type.generated_filename(), 
                    input_in_selectors=self.input_in_selectors,
                    tool=tool,
                    inputs_dict=self.tool.inputs_map()
                )
                if input_type.optional:
                    val = f'{i.id()} && {i.id()} != "None" ? {i.id()} : {val}'
                code = f'def {i.id()} = {val}'
                lines.append(code)

        return lines

    def gen_input_vars(self) -> list[str]:
        """
        Generate Groovy code for input variables definitions inside the Nextflow script section.
        This is where we apply prefix or preprocessiong if necessary.
        """
        lines: list[str] = []
        for a in self.inputs:
            arg_name = ""
            if isinstance(a, ToolInput):
                arg_name = a.id()
            elif isinstance(a, ToolArgument):
                continue
            else:
                raise Exception("unknown input type")

            arg = self.gen_input_var_definition_old(a, arg_name)
            code = f'def {arg_name}WithPrefix = {arg}'  # TODO hack this shit up
            lines.append(code)

        return lines

    def gen_input_var_definition(self, inp: ToolInput, arg_name: str) -> str:
        # def software   = getSoftwareName(task.process)
        # def prefix     = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
        # def read_group = meta.read_group ? "-R ${meta.read_group}" : ""
        # def algoType = params.algoType ? "-a $params.algoType" : ""
        raise NotImplementedError



    def gen_input_var_definition_old(self, inp: ToolInput, arg_name: str) -> str:
        """
        Generate Groovy code to represent the values of input variable definitions
        """
        if isinstance(inp.input_type, Array):
            if (
                isinstance(inp.input_type.subtype(), File)
                and inp.input_type.subtype().has_secondary_files()
            ):
                arg_value = f"get_primary_files({arg_name}).join(' ')"
            else:
                arg_value = f"{arg_name}.join(' ')"

        elif (
            isinstance(inp.input_type, (File)) and inp.input_type.has_secondary_files()
        ):
            arg_value = f"{arg_name}[0]"
        else:
            arg_value = arg_name

        prefix = ""
        if inp.prefix is not None:
            prefix = inp.prefix

        space = ""
        if inp.separate_value_from_prefix is not False:
            space = " "

        prefix_applies_to_all_elements = "False"
        if inp.prefix_applies_to_all_elements is True:
            prefix_applies_to_all_elements = "True"


        # TODO hack this shit up
        if isinstance(inp.input_type, Boolean):
            arg = f"boolean_flag({arg_value}, '{prefix}{space}')"
        else:
            if inp.input_type.optional:
                arg = f"optional({arg_value}, '{prefix}{space}', '{prefix_applies_to_all_elements}')"
            elif inp.prefix:
                #arg = f"apply_prefix({arg_value}, '{prefix}{space}', '{prefix_applies_to_all_elements}')"
                arg = f"'{prefix}{space}' + {arg_value}"
            else:
                arg = arg_value

        return arg





"""
stuff which can be one-liners 
    -A one two three
    -B=seven,eight,nine

stuff which requires function


"""








"""
PREFIX
    -A one two three            return -A ' '.join(items)
    -B=seven,eight,nine         return -B=','.join(items)
    -C=four -C=five -C=six      return ' '.join([f'-C={val}' for val in values])

def apply_prefix(var, prefix, prefix_applies_to_all_elements) {
  
  if (prefix_applies_to_all_elements == 'True') {
      def l = var.split(' ')
    
      prefix = prefix.toString() 
    
      return prefix + l.join(' ' + prefix)
  }
  else {
      return prefix.toString() + var.toString()
  }

}
"""
