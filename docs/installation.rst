============
Installation
============

The package will be available on PyPi. As for now, you can install it by cloning the repo using::
    
    git clone https://github.com/gentangle/pyge

Create a virtual environment and activate it::

    python -m venv env_pyge
    source env_pyge/bin/activate

Navigate inside the folder `pyge` and run::

    python -m pip install .


Recently a convinient tool has been developed to be a similar software to the `cargo` package menager used in Rust, called `uv <https://github.com/astral-sh/uv>`_. Using this tool, the above commands read::

    uv venv env_pyge
    source env_pyge/bin/activate

And in the `pyge` folder::

    uv pip install .

Its convinience it is in its efficiency.