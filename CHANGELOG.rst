
Changelog
=========

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
