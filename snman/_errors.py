class OptionNotImplemented(ValueError):  # pragma: no cover
    """
    Exception for cases when the chosen option is not implemented
    """

    def __init__(self, *args, **kwargs):
        """Create exception"""
        Exception.__init__(self, *args, **kwargs)
