import collections
from sqlalchemy import text, select
from sqlalchemy.dialects.postgresql import insert
import psycopg2.extras
from omnipath_metabo import _log

from . import _structure
from ._base import Base
from ._connection import Connection
from ..misc._tee import Tee

_log("Hello from main module")

def create(con):

    Base.metadata.create_all(con.engine)


TABLES = {'structures', 'identifiers', 'resources'}


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

        loader = Loader(resource, self.con)
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

    def wipe(self) -> None:

        for tbl in TABLES:

            query = text(f'DROP TABLE {tbl} CASCADE')

            self.con.session.execute(query)
            self.con.session.commit()

    def load_all(self):
        for db in _structure.__dict__.values():
            if getattr(db,'scheme',None) is _structure.Structure:
                h = db
                self.load(h)



class Loader():

    def __init__(self, resource, con):
        self._set_resource(resource)
        self.session = con.session
        self.con = con

    def _set_resource(self, resource):
        if type(resource) is type:
            resource = resource()
        self.scheme = resource.scheme
        self.resource = resource

    def load(self):

        self.create()

        insert_resource = insert(_structure.Resource).values(
            name = self.resource.name
        ).returning(_structure.Resource.id)
        insert_resource = insert_resource.on_conflict_do_update(
            index_elements=['name'],
            set_ = {
                'name':insert_resource.excluded.name
            })
        resid = self.session.execute(insert_resource).fetchall()

        _log(f'loading resource {self.resource.name}')

        raw_con = self.con.engine.raw_connection()

        with raw_con.cursor() as cursor:

            query = """
                INSERT INTO structures (name, smiles) VALUES %s ON CONFLICT (smiles) DO nothing;
                """
            _log("loading insert statments for structures table")
            cached_resource = Tee(self.resource, struct = lambda x: x)
            psycopg2.extras.execute_values(cursor, query, cached_resource, page_size = 1000)


        raw_con.commit()
        _log("structures have been inserted, creating mol column")
        return_ids = text("SELECT id, smiles FROM structures")
        inserted_str = set(s[1] for s in cached_resource.struct)
        strids = {
            smiles : id
            for id, smiles in self.session.execute(return_ids)
            if smiles in inserted_str
        }


        #vself.update_mol_column()

        resource_key = resid[0][0]
        _log('resource ids collected.')

        insert_ids = (
            (name, strids[smiles], resource_key)
            for name, smiles in zip(*cached_resource.struct)
        )
        _log('inserting identifiers.')
        with raw_con.cursor() as cursor:
            query = """
                    INSERT INTO identifiers (identifier, structure_id, resource_id) VALUES %s
                    """
            psycopg2.extras.execute_values(cursor, query, insert_ids, page_size = 1000)


        raw_con.commit()
        _log('identifiers inserted.')
        #self.indexer()
        _log(f'{self.resource.name} loaded')

    def update_mol_column(self):
        query = text("update structures set mol = mol_from_smiles(smiles::cstring) where mol is null")
        self.session.execute(query)
        self.session.commit()
        _log('Finished creating mol column.')

    def indexer(self):
        """
        Creates a index of the mol structures using gist. Allows for substructure searches of the molecules.
        """
        query = text("create index molidx on structures using gist(mol)")
        self.session.execute(query)
        self.session.commit()

    def create(self) -> None:

        create(self.con)



