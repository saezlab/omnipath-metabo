from . import _structure
from ._base import Base

def create(con):

    Base.metadata.create_all(con.engine)

class Loader():
    #accept scheme and resource. 
    def __init__(self, scheme, resource, session):
        self.scheme = scheme
        self.resource = resource
        self.session = session


    def load(self):
        for row in self.resource:
            insert_statement = self.scheme.insert().values(
                smiles=row['smiles'],
                accession=row['accession'],
                inchi=row['inchi']
                )
            self.session.execute(insert_statement)
