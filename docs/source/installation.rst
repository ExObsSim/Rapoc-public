=======================
Installation & updates
=======================

Installing from Pypi
--------------------
`RAPOC` can be installed from the Pypi repository with the following script::

    pip install rapoc

Installing from git
-------------------
`RAPOC` may also be cloned from the main git repository::

    git clone https://github.com/ExObsSim/Rapoc-public.git

The next step is to move into the `RAPOC` folder::

    cd /your_path/Rapoc

Then::

    pip install .

To check if one has the correct setup::

    python -c "import rapoc"


Uninstall Rapoc
-------------------

`RAPOC` is installed in your system as a standard python package:
to uninstall it completely::

    pip uninstall rapoc


Update Rapoc
---------------
Update from Pypi
+++++++++++++++++++
To update `RAPOC` to the latest stable version::

    pip install rapoc -upgrade


Update from source
+++++++++++++++++++

One can download or pull a newer version of `RAPOC` over the old one, replacing all modified data.

Then one has to place themselves inside the installation directory within the console::

    cd /your_path/Rapoc

`RAPOC` can now be updated::

    pip install . --upgrade

or simply::

    pip install .


