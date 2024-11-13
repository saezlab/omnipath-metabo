from pypath_common import data
import functools as ft

load = ft.partial(data.load, module = 'omnipath_metabo')