import collections
from sqlalchemy import text, select
from sqlalchemy.dialects.postgresql import insert
import psycopg2.extras
from omnipath_metabo import _log
import itertools

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

    def test_load(self, resource) -> None:

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

        self.chem_properties_populate()
        self.update_mol_column()


# Change the method for updating is_polymer to use InchI.
    def chem_properties_populate(self):
        _log('Querying mol structures', level=-1)
        query = text(r"""
            UPDATE identifiers
            SET
                is_polymer = (structures.inchi ~ '.*\)n.*'),
                has_conformation = (structures.smiles ~ '.*[\\/]=.*'),
                has_stereo = (structures.smiles ~ '.*@.*'),
                complete_formula = NOT (structures.smiles ~ '.*\*.*')
            FROM structures
            WHERE identifiers.structure_id = structures.id;
        """)
        self.con.session.execute(query)
        self.con.session.commit()
        _log('Finished updating columns', level=-1)

    def update_mol_column(self):

        _log("Creating mol column", level =-1)
        q1 = text("update structures set mol = mol_from_smiles(smiles::cstring) where mol is null")
        self.con.session.execute(q1)
        self.con.session.commit()
        _log('Finished creating mol column.', level  =-1)


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
        self.session.commit()

        _log(f'loading resource {self.resource.name}', level = -1)

        raw_con = self.con.engine.raw_connection()

        with raw_con.cursor() as cursor:

            query = """
                INSERT INTO structures (name, smiles, inchi) VALUES %s 
                ON CONFLICT (smiles) 
                DO NOTHING; 
                """
            _log("loading insert statments for structures table")
            cached_resource = Tee(
                self.resource,
                struct = lambda x: x,
                yld = lambda x: x['structure'],
            )
            psycopg2.extras.execute_values(cursor, query, cached_resource, page_size = 1000)


        raw_con.commit()
        _log("structures have been inserted")

        return_ids = text("SELECT id, smiles FROM structures")
        inserted_str = {
            s['structure'][1]
            for s in cached_resource.cached['struct']
        }
        strids = {
            smiles : id
            for id, smiles in self.session.execute(return_ids)
            if smiles in inserted_str
        }


        #vself.update_mol_column()

        resource_key = resid[0][0]
        _log('resource ids collected.')

        insert_ids = (
            (name, strids[smiles], resource_key, True, resource_key)
            for name, smiles, _ in (
                r['structure'] for r in cached_resource.cached['struct']
            )
        )

        _log('inserting identifiers.')
        with raw_con.cursor() as cursor:
            query = """
                    INSERT INTO identifiers (identifier, structure_id, resource_id, authoritative, id_type) VALUES %s
                    ON CONFLICT DO NOTHING
                    """
            psycopg2.extras.execute_values(cursor, query, insert_ids, page_size = 100)


        raw_con.commit()
        _log('identifiers inserted.')

        return_ids = text(
            "SELECT id, structure_id, identifier, resource_id, authoritative, id_type "
            f"FROM identifiers WHERE resource_id = {resource_key}"
        )
        identifier_ids = {
            tuple(rest): id
            for id, *rest in self.session.execute(return_ids)
        }

        property_records = (
            (
                identifier_ids[
                    (
                        strids[record['structure'][1]],
                        record['structure'][0],
                        resource_key,
                        True,
                        resource_key
                     )
                ],
            ) + record['properties']

            for record in cached_resource.cached['struct']
        )

        with raw_con.cursor() as cursor:
            query = """
                    INSERT INTO properties (identifier_id, mw, monoiso_mass, charge, formula) VALUES %s
                    """
            psycopg2.extras.execute_values(cursor, query, property_records, page_size = 1000)


        raw_con.commit()
        _log('properties inserted.')


        #self.indexer()


        _log(f'Finished loading {self.resource.name}.', level = -1)




    def indexer(self):
        """
        Creates a index of the mol structures using gist. Allows for substructure searches of the molecules.
        """
        query = text("create index molidx on structures using gist(mol)")
        self.session.execute(query)
        self.session.commit()

    def create(self) -> None:

        create(self.con)

