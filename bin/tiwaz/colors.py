import matplotlib

import numpy as np


def named_rgb_colors(key, value=None):
    colors = {
        # "blue": (14, 90, 211),
        "blue": (53, 186, 252),
        "green": (23, 114, 6),
        "red": (114, 6, 6),
        "orange": (211, 126, 0),
        "black": (255, 255, 255),
    }
    colors = {k: tuple(c[k] / 256 for k in range(3)) for k, c in colors.items()}

    if value is None:
        return colors[key]

    rgb = colors[key]
    hsv = matplotlib.colors.rgb_to_hsv(rgb)
    hsv[2] = value
    rgb = matplotlib.colors.hsv_to_rgb(hsv)

    return tuple(rgb)


def graded_colors(color_name, n_shades):
    values = np.linspace(0.2, 0.8, n_shades)
    return [named_rgb_colors(color_name, v) for v in values]
