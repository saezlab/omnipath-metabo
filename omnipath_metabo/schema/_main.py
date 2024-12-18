import collections
from sqlalchemy import text, select
from sqlalchemy.dialects.postgresql import insert
import psycopg2.extras
from omnipath_metabo import _log
import itertools
from pypath_common import _misc

from . import _structure
from ._base import Base
from ._connection import Connection
from ..misc._tee import Tee

_log("Hello from main module")

def create(con):

    Base.metadata.create_all(con.engine)


TABLES = {'structures', 'identifiers', 'resources', 'properties'}


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


    def load(self, resource: str, limit: int = None) -> None:

        loader = Loader(resource, self.con)
        loader.load(limit = limit)

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

        for db in list(_structure.__dict__.values()):
            if getattr(db,'scheme',None) is _structure.Structure:
                h = db
                self.load(h)

        self.chem_properties_populate()
        #self.update_mol_column()


# Change the method for updating is_polymer to use InchI.
    def chem_properties_populate(self):
        _log('Querying mol structures', level=-1)
        query = text(r"""
            UPDATE structures
            SET
                is_polymer = (inchi ~ '.*\\p.*'),
                has_conformation = (smiles ~ '.*[\\/]=.*'),
                has_stereo = (smiles ~ '.*@.*'),
                complete_formula = NOT (smiles ~ '.*\*.*')
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
        self._resource_ids = {}

    def _set_resource(self, resource: str | type):
        if type(resource) is type:
            resource = resource()
        self.scheme = resource.scheme
        self.resource = resource

    def load(self, limit: int = None):

        self.create()

        resource_key = self._resource_id(self.resource.name)


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
                itertools.islice(self.resource, limit),
                struct = lambda x: x,
                yld = lambda x: x['structure'],
            )
            psycopg2.extras.execute_values(cursor, query, cached_resource, page_size = 1000)

        type(cached_resource)

        raw_con.commit()
        _log("structures have been inserted", level = -1)

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
        _log('resource ids collected.', level = -1)

        insert_ids = (
            (
                _id,
                strids[smiles],
                resource_key,
                id_type == self.resource.name,
                self._resource_id(self.resource.id_types.get(id_type, id_type))
            )
            for i, r in enumerate(cached_resource.cached['struct'])
            for name, smiles, _ in [r['structure']]
            for _id, id_type in itertools.chain(
                [(name, self.resource.name)],
                r['identifiers']
            )
        )

        _log('inserting identifiers.', level =-1)
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


    def _resource_id(self, name: str | list[str]) -> int:

        name = _misc.to_set(name)

        if new_names := sorted(name - set(self._resource_ids.keys())):

            insert_resource = (
                insert(_structure.Resource).
                values([{'name': l} for l in new_names]).
                returning(_structure.Resource.id)
            )

            insert_resource = insert_resource.on_conflict_do_update(
                    index_elements=['name'],
                    set_ = {
                        'name':insert_resource.excluded.name
                    }
                )

            resids = self.session.execute(insert_resource).fetchall()
            self.session.commit()

            self._resource_ids.update({
                label: _id[0]
                for label, _id in zip(new_names, resids)
            })

        return (
            [self._resource_ids[n] for n in name]
                if len(name) > 1 else
            self._resource_ids[_misc.first(name)]
        )


    def indexer(self):
        """
        Creates a index of the mol structures using gist. Allows for substructure searches of the molecules.
        """
        query = text("create index molidx on structures using gist(mol)")
        self.session.execute(query)
        self.session.commit()

    def create(self) -> None:

        create(self.con)


    def _id_type(self, id_idx: int) -> int:

        self.resource.name
        self.resource.id_fields[id_idx]


