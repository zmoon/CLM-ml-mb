def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark tests as slow (deselect with `-m 'not slow'`)")
