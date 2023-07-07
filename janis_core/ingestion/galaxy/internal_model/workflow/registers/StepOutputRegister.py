


from ..step.outputs import StepOutput


class StepOutputRegister:
    def __init__(self):
        self.register: list[StepOutput] = []

    def add(self, step_output: StepOutput) -> None:
        self.register.append(step_output)

    def list(self) -> list[StepOutput]:
        return self.register
        

