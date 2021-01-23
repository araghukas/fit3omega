from setuptools import setup, find_packages, Extension
from fit3omega import __version__

C_module = Extension('intglib', sources=['./integrate/intglib.c',
                                         './integrate/intgmodule.c'])

setup(
    name="fit3omega",
    version=__version__,
    description="Data analyzer for 3ω thermal conductivity measurements",
    author="Ara Ghukasyan",
    author_email="ghukasa@mcmaster.ca",
    license="MIT License",
    packages=find_packages(),
    install_requires=['pyyaml', 'pandas', 'numpy', 'matplotlib'],
    ext_modules=[C_module]
)
