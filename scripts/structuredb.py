from omnipath_metabo.schema import _main
from omnipath_metabo.schema import _structure

d = {
    'user': 'metabo',
    'password': 'metabo123',
    'host': 'localhost',
    'port': '5432',
    'database': 'metabo',
}
db = _main.Database(d)
db.wipe()
ramp = _structure.Ramp()
# db.load_all()
db.load(ramp)
