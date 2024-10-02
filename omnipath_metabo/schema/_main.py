from sqlalchemy import insert

from . import _structure
from ._base import Base
from ._connection import Connection

def create(con):

    Base.metadata.create_all(con.engine)


class Database:

    def __init__(self, con):

        self.con = con
        self.connect()
        self.create()


    def connect(self, reconnect: bool = False) -> None:

        if reconnect or not isinstance(self.con, Connection):

            self.con = Connection(**self.con)
            self.con.connect()


    def create(self) -> None:

        create(self.con)


    def load(self, resource) -> None:

        loader = Loader(resource, self.con.session)
        loader.load()


    def __del__(self):

        if hasattr(self, 'con'):

            del self.con


class Loader():
    def __init__(self, resource, session):
        self.scheme = resource.scheme
        self.resource = resource
        self.session = session


    def load(self):

        for row in self.resource:

            insert_statement = insert(self.scheme).values(
                structure=row[1],
                name=row[0],
                #inchi=row['inchi']
                )

            self.session.execute(insert_statement)

        self.session.commit()
