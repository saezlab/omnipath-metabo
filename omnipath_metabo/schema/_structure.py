from sqlalchemy import (
    Boolean,
    create_engine,
    Column,
    ForeignKey,
    Integer,
    String,
    types,
    Numeric,
    Index,
    UniqueConstraint,
    func
)
from sqlalchemy.orm import relationship
from ._base import Base
from pypath.inputs import hmdb, swisslipids, lipidmaps, ramp
from .. import data as _data

class MolType(types.UserDefinedType):
    cache_ok = True

    def get_col_spec(self, **kw):
        return "MOL"


class Structure(Base):
    __tablename__ = 'structures'
    id = Column(Integer, primary_key = True)
    name = Column(String)
    smiles = Column(String, unique= True)
    inchi = Column(String(4096))
    mol = Column(MolType)
    identifier = relationship(
        'Identifier',
        backref='structure',
        foreign_keys='Identifier.structure_id'
    )

class Identifier(Base):
    __tablename__ = 'identifiers'
    id = Column(Integer, primary_key=True)
    identifier = Column(String)
    structure_id = Column(
        Integer,
        ForeignKey('structures.id'),
        nullable = False,
    )
    resource_id = Column(
        Integer,
        ForeignKey('resources.id'),
        nullable = False
    )
    authoritative = Column(Boolean)
    id_type = Column(
        Integer,
        ForeignKey('resources.id'),
        nullable = False)
    is_polymer = Column(Boolean)
    has_stereo = Column(Boolean)
    has_conformation = Column(Boolean)
    complete_formula = Column(Boolean)
    __table_args__ = (
        Index(
            'identifers_unique',
            'identifier',
            'structure_id',
            'resource_id',
            'id_type',
            unique = True,
        ),
    )
    properties = relationship(
        'Properties',
        backref='properties',
        foreign_keys='Properties.identifier_id'
    )



class Resource(Base):
    __tablename__ = 'resources'
    id = Column(Integer, primary_key=True)
    name = Column(String, unique = True, nullable = False)
    identifier = relationship(
        'Identifier',
        backref='resource',
        foreign_keys='Identifier.resource_id'
    )

class Properties(Base):
    __tablename__ = 'properties'
    id = Column(Integer, primary_key=True)
    identifier_id = Column(Integer, ForeignKey('identifiers.id'), nullable = False)
    mw = Column(Numeric)
    monoiso_mass = Column(Numeric)
    charge = Column(Numeric)
    formula = Column(String)

class Mappings(Base):
    __tablename__ = 'mappings'
    id = Column(Integer, primary_key= True)
    structure_id = Column(Integer, ForeignKey('structures.id'), nullable = False)
    resource_id = Column(String)
    chebi_id = Column(String)
    pubchem_id = Column(String)



class ResourceBase:
    def __init__(self):
        self.id_types = self._resource_label()

    @classmethod
    def _resource_label(cls) -> str:

        return [
            _id_types().get(cls.name, {}).get(label, label)
            for label in cls.id_fields
        ]



class Hmdb(ResourceBase):
    scheme = Structure
    name = 'HMDB'
    id_fields = sorted(f'{s}_id' for s in hmdb.ID_FIELDS)

    def __iter__(self):
        for met in hmdb.metabolites_processed(
            'accession',
            'smiles',
            'average_molecular_weight',
            'monisotopic_molecular_weight',
            'chemical_formula',
            'inchi',
            *self.id_fields,
            ):

            yield {
                'structure':(met[0][0], met[0][1], met[0][5]),
                'properties': (met[0][2], met[0][3], None, met[0][4]),
                'identifiers': met[0][6:]
            }


class SwissLipids():
    scheme = Structure
    name = 'SwissLipids'
    def __iter__(self):
        for met in swisslipids.swisslipids_lipids():
            if met['Mass (pH7.3)'] == '':
                mass = 0
            else:
                mass = met['Mass (pH7.3)']
            if met['Charge (pH7.3)'] == '':
                charge = 0
            else:
                charge = met['Charge (pH7.3)']

            yield {
                'structure': (
                    met['Lipid ID'],
                    met['SMILES (pH7.3)'],
                    met['InChI (pH7.3)'],
                ),
                'properties':(mass, 0, charge, met['Formula (pH7.3)'])
                }


class LipidMaps():
    scheme = Structure
    name = 'LipidMaps'
    def __iter__(self):
        sdf = lipidmaps.lmsd_sdf()

        for met in sdf:

            if smiles := met['name'].get('SMILES', None):

                if met['id'] == 'LMFA08040060':

                    continue

                yield {
                    'structure':(met['id'], smiles, met['name']['INCHI']),
                    'properties':(met['annot'].get('EXACT_MASS', None), None, None, met['name'].get('FORMULA', None) )
                }

class Ramp(ResourceBase):
    scheme = Structure
    name = 'RaMP'
    id_fields = sorted(ramp.ramp_id_types())

    def __iter__(self):

        for row in ramp.ramp_iter('chem_props'):

            yield {
                'structure': (row[0], row[3], row[6]),
                'properties':(row[7], row[8], None, row[10])
                }


def _id_types() -> dict[str, dict]:

    if 'ID_TYPES' not in globals():
        globals()['ID_TYPES'] = _data.load('resource-labels')

    return globals()['ID_TYPES']
