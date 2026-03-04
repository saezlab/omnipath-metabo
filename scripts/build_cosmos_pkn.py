#!/usr/bin/env python

"""
Thin wrapper — forwards to the installed ``cosmos-pkn`` entry point.

Prefer calling the entry point directly::

    cosmos-pkn [options]

or via uv::

    uv run cosmos-pkn [options]

See ``cosmos-pkn --help`` for all options.
"""

from omnipath_metabo.datasets.cosmos._cli import main

if __name__ == '__main__':
    main()
