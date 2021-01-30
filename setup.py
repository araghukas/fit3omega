from setuptools import setup, find_packages, Extension
from fit3omega import __version__

C_module = Extension('intg', sources=['./integrate/intg.c'])

#
# C_module = Extension('integrate',
#                      sources=["./integrate/integrate.cpp"],# "./integrate/module.cpp"],
#                      language="c++",
#                      extra_compile_args=["-std=c++11", "-stdlib=libc++"],
#                      extra_link_args=["-stdlib=libc++"])

setup(
    name="fit3omega",
    version=__version__,
    description="Data analyzer for 3ω thermal conductivity measurements",
    author="Ara Ghukasyan",
    author_email="ghukasa@mcmaster.ca",
    license="MIT License",
    packages=find_packages(),
    package_data={"fit3omega": ["integrate/*"]},
    install_requires=['pyyaml', 'pandas', 'numpy', 'matplotlib'],
    ext_modules=[C_module]
)
