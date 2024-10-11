from sqlalchemy import insert, text, select

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

    def substructure_search(self, substructure):

        query = f"select name, mol from structures where mol @>'{substructure}'"
        result = self.con.session.execute(text(query))
        for row in result:
            print(f"{row[0]}, {row[1]}")
        return result

    def __del__(self):

        if hasattr(self, 'con'):

            del self.con


class Loader():
    def __init__(self, resource, session):
        self.scheme = resource.scheme
        self.resource = resource
        self.session = session


    def load(self):

        for i, row in enumerate(self.resource):

            insert_statement = insert(self.scheme).values(
                smiles=row[1],
                name=row[0],
                
                )

            self.session.execute(insert_statement)
            if i > 1000:
                break
        self.session.commit()
        self.update_mol_column()
        self.indexer()

    def update_mol_column(self):
        query = text("update structures set mol = mol_from_smiles(smiles::cstring) where mol is null")
        self.session.execute(query)
        self.session.commit()

    def indexer(self):
        """
        Creates a index of the mol structures using gist. Allows for substructure searches of the molecules.
        """
        query = text("create index molidx on structures using gist(mol)")
        self.session.execute(query)
        self.session.commit()



        
