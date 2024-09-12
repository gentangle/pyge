__version__ = "0.8.9"

import sys

from loguru import logger

logger.remove(0)
logger.add(sys.stderr, level="INFO")
