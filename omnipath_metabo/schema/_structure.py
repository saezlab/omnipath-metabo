from sqlalchemy import create_engine, Column, Integer, String
from ._base import Base


class Structure(Base):
    __tablename__ = 'structures'
    id = Column(Integer, primary_key=True)
    name = Column(String)
    structure = Column('mol_structure', String)
