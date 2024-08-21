import argparse

class ListAppend(argparse.Action):
    """Custom action for appending values to a list.

    This action ensures that the target attribute in the namespace is initialized
    as a list before appending values.

    Args:
        parser: The ArgumentParser instance.
        namespace: The namespace object where the parsed arguments are stored.
        values: A list of values to be appended.
        option_string: The string representation of the option (e.g., '-o', '--option').
    """
    def __call__(self, parser, namespace, values, option_string=None):
        if values:
            if getattr(namespace, self.dest) is self.default:
                setattr(namespace, self.dest, [])
            items = getattr(namespace, self.dest, None)
            setattr(namespace, self.dest, sorted(set(items + values)))