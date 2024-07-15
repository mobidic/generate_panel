import setuptools
import os

def get_version():
    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'VERSION')) as version_file:
        version = version_file.read().strip()
    return version

setuptools.setup(
    version=get_version()
)