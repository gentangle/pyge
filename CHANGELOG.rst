
Changelog
=========

v0.8.9 (2024-09-12)
------------------------------------------------------------

v0.8.8 (2024-03-13)
------------------------------------------------------------

* Switching from `logging` to `loguru` for logging
* Faster frame section when using `trajectory`
* Now `trajectory` can handle more `mask` options to select which frame use to compute GE: a list of frames, a numpy bool with selected frames, a dictionary with `start`, `stop` and `step` keys
* `trajectory` now is consistent with `singlechain`: both automatically computes the contact map for the given topology 

v0.8.7 (2024-03-13)
------------------------------------------------------------

* Doc update: installation section revamped.
* New minimal dependency requirements
* Simplified linting dependency: using `ruff` for linting (see `tox.ini`). This is the raccomanded choice to format code too.

v0.8.6 (2023-05-29)
------------------------------------------------------------

v0.8.5 (2023-05-26)
------------------------------------------------------------

* update documetation
* specify with lower bounds libraries versions

v0.8.4 (2023-05-26)
------------------------------------------------------------

v0.8.3 (2023-05-05)
------------------------------------------------------------

* `pyproject.toml` added for pre-build specifications (Cython)

v0.8.2 (2023-03-29)
------------------------------------------------------------

v0.8.1 (2023-03-29)
------------------------------------------------------------

v0.8.0 (2023-03-28)
------------------------------------------------------------

* PDB parser now has an option to avoid `HETATM` records

v0.7.0 (2023-03-27)
------------------------------------------------------------

* adapt native contact map function to address those sequences that presents gaps.

v0.6.1 (2023-03-21)
------------------------------------------------------------

v0.6.0 (2023-03-20)
------------------------------------------------------------

* Implement three new objects to store GE calculations: `GE`, `GETermini` and `GEChain`. `GE` collects the loop, the thread and the GE value of a loop calculation. `GETermini` object contains the `GE`` for the N- and C-threads. `GEChain` collects the output of `pyge.singlechain.singlechain`

v0.5.0 (2023-03-20)
------------------------------------------------------------

* Integrate the computation of GE weighted (see publication) into the main function in `gent`

v0.4.0 (2023-03-20)
------------------------------------------------------------

v0.3.0 (2023-03-20)
------------------------------------------------------------

v0.2.0 (2023-03-20)
------------------------------------------------------------

v0.1.1 (2023-03-20)
------------------------------------------------------------

* Implementation of a new `singlechain` function. Goal: provide a function to return GE values directly from a PDB file
* Remove `qentangled` module because it is not used anymore

v0.1.0 (2022-12-15)
------------------------------------------------------------

* First release
