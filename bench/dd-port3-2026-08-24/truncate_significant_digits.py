#!/usr/bin/env python3
"""Truncate every numeric token in an SDPA sparse file to N significant digits.

The point is to vary TOKEN LENGTH while leaving the value alone at the target precision. dd_real
carries about 32 significant digits, so N well above 32 (and far below QD's ~309-digit converter
limit) changes the text without changing what a correct reader produces.
"""
import sys, re
from decimal import Decimal, getcontext

n = int(sys.argv[2])
getcontext().prec = n
num = re.compile(r'[-+]?\d+\.\d+(?:[eE][-+]?\d+)?')

def shrink(mo):
    d = Decimal(mo.group(0))
    # +0 applies the context precision, i.e. rounds to n significant digits.
    return format(d + Decimal(0), 'f')

with open(sys.argv[1]) as fh:
    for line in fh:
        sys.stdout.write(num.sub(shrink, line))
