from sqlalchemy import create_engine, Column, Integer, String
from ._base import Base
from pypath.inputs import hmdb


class Structure(Base):
    __tablename__ = 'structures'
    id = Column(Integer, primary_key=True)
    name = Column(String)
    structure = Column('mol_structure', String)


class Hmdb():
    scheme = Structure

    def __iter__(self):

        for met in hmdb.metabolites_processed('name','smiles'):
            yield met[0]

class SwissLipids():
    scheme = Structure 

    def __iter__(self):
        for met in swisslipids.metabolites('name','smiles'):
            yield met