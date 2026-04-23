# Shared list utilities.
# Andrew Benson (ported to Python 2026)


def as_array(item):
    """Normalize None/scalar/list → list (port of List::ExtraUtils::as_array)."""
    if item is None:
        return []
    if isinstance(item, list):
        return item
    return [item]
