
import os
import yaml
from contextlib import closing

from sqlalchemy import create_engine

class Connection:
    def __init__(
            self,
            param: str | dict | None = None,
            **kwargs
        ):

        self._param = param or kwargs
        self._parse_param()


    def _parse_param(self) -> None:

        self._from_file()


    def _from_file(self) -> None:

        if isinstance(self._param, str) and os.path.exists(self._param):

            with closing(open(self._param, 'r')) as fp:

                self._param = yaml.load(fp, Loader = yaml.FullLoader)

    def _uri(self) -> str:

        return f'postgresql://'
