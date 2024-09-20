from omnipath_metabo.schema import _connection

def test_param_yaml(dummy_connection_config, dummy_con_param):

    con = _connection.Connection(dummy_connection_config)

    assert con._param == dummy_con_param

def test_param_dict(dummy_con_param):

    con = _connection.Connection(dummy_con_param)

    assert con._param == dummy_con_param