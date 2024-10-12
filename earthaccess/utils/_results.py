from functools import reduce
from typing import Dict, Union
from ipyleaflet import Map, Marker, Rectangle
from ipywidgets import Output


def get_nested_value(data, keys, default=None) -> Union[Dict, None]:
    """Safely access a nested value in a dictionary."""
    try:
        return reduce(lambda d, key: d[key], keys, data)
    except (KeyError, IndexError, TypeError):
        return default
    
def convert_to_mb(size: float, unit: str) -> float:
        if unit == "MB":
            return size
        elif unit == "GB":
            return size * 1024
        elif unit == "KB":
            return size / 1024
        else:
            return size / (1024 * 1024)
