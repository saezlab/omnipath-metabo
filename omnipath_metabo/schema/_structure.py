from sqlalchemy import create_engine, Column, Integer, String, types
from ._base import Base
from pypath.inputs import hmdb, swisslipids

class MolType(types.UserDefinedType):
    cache_ok = True

    def get_col_spec(self, **kw):
        return "MOL"


class Structure(Base):
    __tablename__ = 'structures'
    id = Column(Integer, primary_key=True)
    name = Column(String)
    smiles = Column(String)
    mol = Column(MolType)


class Hmdb():
    scheme = Structure

    def __iter__(self):

        for met in hmdb.metabolites_processed('name','smiles'):
            yield met[0]

class SwissLipids():
    scheme = Structure

    def __iter__(self):
        for met in swisslipids.swisslipids_lipids():
            yield met['Lipid ID'], met['SMILES (pH7.3)']

#Ramp and LipidMaps
