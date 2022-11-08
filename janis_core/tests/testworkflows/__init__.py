

from .core_features import (
    
    # basics
    BasicIOTestWF,
    StepInputsTestWF,
    StepInputsWFInputTestWF,
    StepInputsStaticTestWF,
    StepInputsPartialStaticTestWF,
    StepInputsMinimalTestWF,
    StepConnectionsTestWF,

    # arrays
    ArrayIOTestWF,
    ArrayStepInputsTestWF,
    ArrayStepConnectionsTestWF,

    # scatter
    BasicScatterTestWF,
    ChainedScatterTestWF,
    MultiFieldScatterTestWF,

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

from .assembly import w as AssemblyTestWF