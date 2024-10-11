from sqlalchemy import create_engine, Column, ForeignKey, Integer, String, types
from sqlalchemy.orm import relationship
from ._base import Base
from pypath.inputs import hmdb, swisslipids, lipidmaps

class MolType(types.UserDefinedType):
    cache_ok = True

    def get_col_spec(self, **kw):
        return "MOL"


class Structure(Base):
    __tablename__ = 'structures'
    id = Column(Integer, primary_key=True)
    name = Column(String)
    smiles = Column(String, unique = True)
    mol = Column(MolType)
    identifier = relationship('Identifier', backref='structure')


class Identifier(Base):
    __tablename__ = 'identifiers'
    id = Column(Integer, primary_key=True)
    identifier = Column(String, unique = True)
    structure_id = Column(
        Integer,
        ForeignKey('structures.id'),
        nullable = False,
    )
    resource_id = Column(Integer, ForeignKey('resources.id'), nullable = False)


class Resource(Base):
    __tablename__ = 'resources'
    id = Column(Integer, primary_key=True)
    name = Column(String)
    identifier = relationship('Identifier', backref='resource')


class Hmdb():
    scheme = Structure
    name = 'HMDB'
    def __iter__(self):

        for met in hmdb.metabolites_processed('name','smiles'):
            yield met[0]

class SwissLipids():
    scheme = Structure
    name = 'SwissLipids'
    def __iter__(self):
        for met in swisslipids.swisslipids_lipids():
            yield met['Lipid ID'], met['SMILES (pH7.3)']
"""
class LipidMaps():
    scheme = Structure

    def __iter__(self):
        sdf = lipidmaps.structure.
        for met in sdf.data:
            yield met['LM_ID'], met['SMILES']
"""

#Ramp and LipidMaps
