from logging import getLogger

logger = getLogger(__name__)


def set_log_level(level):
    '''Set the debug level for calour

    You can see the logging levels at:
    https://docs.python.org/3.5/library/logging.html#levels

    Parameters
    ----------
    level : int or str
        10 for debug, 20 for info, 30 for warn, etc.
        It is passing to :func:`logging.Logger.setLevel`

    '''
    logger = getLogger(__package__)
    logger.setLevel(level.upper())
