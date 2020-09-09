import pkg_resources


def get_data():
    return pkg_resources.resource_filename(__name__, 'data/')
