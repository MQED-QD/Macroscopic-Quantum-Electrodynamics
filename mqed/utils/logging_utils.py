import logging
from loguru import logger
import sys

# This InterceptHandler is a standard pattern for bridging the two logging systems.
class InterceptHandler(logging.Handler):
    """
    Takes messages from the standard logging module and redirects them to Loguru.
    """
    def emit(self, record):
        try:
            level = logger.level(record.levelname).name
        except ValueError:
            level = record.levelno

        frame, depth = logging.currentframe(), 2
        while frame.f_code.co_filename == logging.__file__:
            frame = frame.f_back
            depth += 1

        logger.opt(depth=depth, exception=record.exc_info).log(level, record.getMessage())

def setup_loguru():
    """
    Configures Loguru to be the single, authoritative logger.
    It removes default handlers, adds a new one that respects Hydra's
    console logging level, and intercepts standard logging messages.
    """
    # Remove the default loguru handler
    logger.remove()
    
    # Get the logging level from Hydra's config (or default to INFO)
    log_level = logging.getLogger().level or logging.INFO
    
    # Add a new sink to stderr that respects the configured level
    logger.add(sys.stderr, level=log_level)
    
    # Add the intercept handler to capture standard logging messages
    logging.basicConfig(handlers=[InterceptHandler()], level=0, force=True)