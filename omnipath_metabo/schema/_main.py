from sqlalchemy.ext.declarative import declarative_base


__all__ = ['Base']

Base = declarative_base()


def create(con):

    Base.metadata.create_all(con.engine)
