import os
import pytest
import yaml

@pytest.fixture
def testtmp(tmpdir_factory):

    return tmpdir_factory.mktemp('omnipath-metabo-tests')

@pytest.fixture
def dummy_con_param():

    return {
        'user': 'testuser',
        'password': 'testpw',
        'host': 'myhost',
        'port': 5432,
        'database': 'testdb',
    }



@pytest.fixture
def dummy_connection_config(testtmp, dummy_con_param):

    path = os.path.join(testtmp, 'config.yaml')

    with open(path, 'w') as fp:

        fp.write(yaml.dump(dummy_con_param))

    return path
