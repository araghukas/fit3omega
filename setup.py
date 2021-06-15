import numpy as np
from setuptools import setup, find_packages, Extension
from fit3omega import __version__

C_module = Extension("integrate",
                     sources=["./integrate/integrate.c",
                              "./integrate/util.c",
                              "./integrate/borca_tasciuc.c",
                              "./integrate/olson_graham_chen.c"],
                     include_dirs=[np.get_include()])

setup(
    name="fit3omega",
    version=__version__,
    description="Data analyzer for 3ω thermal conductivity measurements",
    author="Ara Ghukasyan",
    author_email="ghukasa@mcmaster.ca",
    license="MIT License",
    packages=find_packages(),
    package_data={"fit3omega": ["integrate/*"]},
    install_requires=['pyyaml', 'pandas', 'numpy', 'matplotlib', 'scipy'],
    ext_modules=[C_module]
)
