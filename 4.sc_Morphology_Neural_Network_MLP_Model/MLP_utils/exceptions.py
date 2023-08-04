# define python user-defined exceptions


class TrainingValidationTestingSplitError(Exception):
    """Exception occurs in the case that the training, validation, and testing splits do not add up to 1.0"""

    pass


class YDataTypeError(Exception):
    """Exception occurs in the case that the y data does not have the correct data type
    (e.g., each model has a different y data type)"""

    pass


class ModelTypeError(Exception):
    """Exception occurs in the case that the model type is not specified correctly"""

    pass


class OptimizationMetricError(Exception):
    """Exception occurs in the case that the optimization metric is not specified correctly"""

    pass


class ModelNameError(Exception):
    """Exception occurs in the case that the model name is not specified correctly"""

    pass
