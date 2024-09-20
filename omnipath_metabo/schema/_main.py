from . import _structure
from ._base import Base

def create(con):

    Base.metadata.create_all(con.engine)
