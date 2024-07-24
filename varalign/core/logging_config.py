# varalign/core/logging_config.py

import logging
import sys


def configure_logging():
    logging.getLogger("varalign").addHandler(logging.NullHandler())
    logging.captureWarnings(True)
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
    )
