Top-level notebooks are run without the library to demonstrate a static version of each of the components.

Graphs are saved in the "Outputs" directory.

Velocity arrays are saved in the "Inputs" directory. The structure of these are as follows:
    - Filetype: .hdf5
    - Filename: varray_rmin-rmax_lenr.hdf5
    - Group: "bulge", "blackhole", "disk", "halo"
    - Dataset: characteristic_value# (ex: n2.7, h3, Mbh12)
Scripts should automatically sort datsets according to this naming convention so that the correct one can be pulled.
