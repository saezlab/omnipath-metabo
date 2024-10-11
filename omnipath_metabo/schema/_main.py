import collections
from sqlalchemy import text, select
from sqlalchemy.dialects.postgresql import insert

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

        insert_resource = insert(_structure.Resource).values(
            name = self.resource.name
        )
        insert_resource = insert_resource.on_conflict_do_nothing(index_elements=['name'])
        self.session.execute(insert_resource)
        ids = collections.defaultdict(set)

        for i, row in enumerate(self.resource):

            insert_statement = insert(self.scheme).values(
                smiles=row[1],
                name=row[0],
            )
            ids[row[1]].add(row[0])

            insert_statement = insert_statement.on_conflict_do_nothing(index_elements = ['smiles'])
            self.session.execute(insert_statement)

            if i > 1000:
                break

        self.session.commit()
        self.update_mol_column()

        select_str_ids = text('SELECT id, smiles FROM structures')
        strids = {
            id[1]: id[0]
            for id in self.session.execute(select_str_ids)
        }
        
        select_res_ids = text('SELECT id, name FROM resources')
        resid= {
            id[1]: id[0]
            for id in self.session.execute(select_res_ids)
        }
        resource_key = resid[self.resource.name]
        insert_ids = insert(_structure.Identifier).values([
            {'identifier':id, 'structure_id': strids[smiles], 'resource_id': resource_key}
            for smiles, _ids in ids.items()
            for id in _ids
        ])
        self.session.execute(insert_ids)
        self.session.commit()
    

        #self.indexer()

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




