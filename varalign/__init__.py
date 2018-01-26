#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import sys

logging.getLogger("varalign").addHandler(logging.NullHandler())
logging.captureWarnings(True)
logging.basicConfig(stream=sys.stderr, level=logging.INFO,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s')


__title__ = 'varalign'
__version__ = '0.0.0'
__license__ = 'MIT'
__authors__ = 'Stuart MacGowan'
