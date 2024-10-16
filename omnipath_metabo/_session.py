import functools as _ft

from pypath_common import log as _read_log
from pypath_common import session as _session

_get_session = _ft.partial(_session, 'omnipath_metabo')
log = _ft.partial(_read_log, 'omnipath_metabo')

session = _get_session()
_log = session._logger.msg
