from sqlalchemy import insert, text

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

    def update_mol_column(self):
        query = text("update structures set mol = mol_from_smiles(smiles::cstring) where mol is null")
        self.session.execute(query)
        self.session.commit()

        


"""
After creating the entries into the database, should create rdkit extension. 
Then a molindexer should be created to allow for efficient substructure searching. 
The psql commands for this are 
CREATE EXTENSION IF NOT EXISTS rdkit ;
-CREATE
CREATE SCHEMA rdk;
-CREATE
select * into mols from (select id,mol_from_smiles(smiles::cstring) m from raw_data) tmp where m is not null;
-SELECT 270010
CREATE INDEX molidx ON mols USING gist(m);
-CREATE INDEX 

One problem is that the 


"""

        
