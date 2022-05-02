# fit3omega
A Python package for fitting voltage data from 3-omega thermal conductivity measurements.
This is a work in progress.

## installation
In the directory that contains `setup.py`:

    pip install .

this automatically should compile and link the C-extension module.

## usage
`fit3omega` is useful once 3-omega measurement and error data has been acquired, and
a sample configuration file prepared (see example folder). The quickest way to run a data fit
is to invoke `fit3omega` from the command line:

    python -m fit3omega sample.txt data.csv

Try it in the `example/` directory.
