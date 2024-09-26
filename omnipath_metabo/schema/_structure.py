from sqlalchemy import create_engine, Column, Integer, String
from ._base import Base
from pypath.inputs import hmdb

class Structure(Base):
    __tablename__ = 'structures'
    id = Column(Integer, primary_key=True)
    name = Column(String)
    structure = Column('mol_structure', String)


class Hmdb(Base):
    scheme = Structure

    def __init__(self):

        pass


    def __iter__(self):

        for met in hmdb.iter_metabolites():
