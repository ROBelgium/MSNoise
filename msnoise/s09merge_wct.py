#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Backward-compatibility shim — replaced by s09compute_wct_dtt.py.
"""
from .s09compute_wct_dtt import main  # noqa: F401


def main(loglevel="INFO", wct_dir=None, output_dir=None):
    from .s09compute_wct_dtt import main as _main
    _main(loglevel=loglevel)


if __name__ == "__main__":
    main()
