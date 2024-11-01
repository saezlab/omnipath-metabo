from sqlalchemy import Boolean, create_engine, Column, ForeignKey, Integer, String, types, Numeric, Index
from sqlalchemy.orm import relationship
from ._base import Base
from pypath.inputs import hmdb, swisslipids, lipidmaps, ramp

class MolType(types.UserDefinedType):
    cache_ok = True

    def get_col_spec(self, **kw):
        return "MOL"


class Structure(Base):
    __tablename__ = 'structures'
    id = Column(Integer, primary_key=True)
    name = Column(String)
    smiles = Column(String, unique= True)
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
    id_type = Column(Integer, ForeignKey('resources.id'), nullable = False)
    __table_args__ = (Index('identifers_unique','identifier', 'structure_id', 'resource_id', 'id_type', unique = True),)



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


class Hmdb():
    scheme = Structure
    name = 'HMDB'
    def __iter__(self):

        for met in hmdb.metabolites_processed(
            'accession',
            'smiles'
            ):
            yield met[0]

class SwissLipids():
    scheme = Structure
    name = 'SwissLipids'
    def __iter__(self):
        for met in swisslipids.swisslipids_lipids():
            yield met['Lipid ID'], met['SMILES (pH7.3)']


class LipidMaps():
    scheme = Structure
    name = 'LipidMaps'
    def __iter__(self):
        sdf = lipidmaps.lmsd_sdf()

        for met in sdf:

            if smiles := met['name'].get('SMILES', None):

                if met['id'] == 'LMFA08040060':

                    continue

                yield met['id'], smiles

class Ramp():
    scheme = Structure
    name = 'RaMP'
    def __iter__(self):

        for row in ramp.ramp_iter('chem_props'):

            yield {'structure': (row[0], row[3]), 'properties':(row[7], row[8], None, row[10])}
