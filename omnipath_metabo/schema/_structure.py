from sqlalchemy import create_engine, Column, Integer, String
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

class Structure(Base):
    __tablename__ = 'structures'
    id = Column(Integer, primary_key=True)
    mol_name = Column(String)
    mol_structure = Column('mol_structure', String)
