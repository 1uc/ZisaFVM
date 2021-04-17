def with_default(primary_value, default_value):
    """
    Return `primary_value` unless it's `None`.

    Examples:
        >>> def foo(bar=None):
        >>>     bar = with_default(bar, 42.0)

        >>> print(foo())
        >>> 42.0
        >>> print(foo("Convenient."))
        >>> Convenient.
    """

    if primary_value is not None:
        return primary_value
    else:
        return default_value
