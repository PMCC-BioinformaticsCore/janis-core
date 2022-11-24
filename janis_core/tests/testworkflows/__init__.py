

from .core_features import (
    
    # basics
    BasicIOTestWF,
    WildcardSelectorOutputTestWF,
    InputSelectorOutputTestWF,
    StepInputsTestWF,
    StepInputsWFInputTestWF,
    StepInputsStaticTestWF,
    StepInputsPartialStaticTestWF,
    StepInputsMinimalTestWF,
    StepConnectionsTestWF,
    DirectivesTestWF,

    # arrays
    ArrayIOTestWF,
    ArrayIOExtrasTestWF,
    ArrayStepInputsTestWF,
    ArrayStepConnectionsTestWF,

    # scatter
    BasicScatterTestWF,
    ChainedScatterTestWF,
    ScatterDotTestWF,
    ScatterCrossTestWF,

    # secondaries
    SecondariesIOTestWF,
    SecondariesConnectionsTestWF,

    # combos
    ScatterSecondariesTestWF,
    ArraySecondariesTestWF

)

from .additional_features import (
    StepInputExpressionTestWF,
    ConditionStepTestWF,
    AliasSelectorTestWF,
    ArraysOfSecondaryFilesOutputsTestWF,
    ForEachTestWF
)

from .codetools import (
    InputsPythonToolTestWF,
    OutputsPythonToolTestWF
)

from .assembly import w as AssemblyTestWF