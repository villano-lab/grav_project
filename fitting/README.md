# Fitting

This directory contains:
- All notebooks in which we curve fit to data for NGC5533
- Outputs of curve fitting attempts
- .hdf5 files with pre-calculated rotation curves to be used for faster fitting routines. This is especially important for the disk calculation, which takes a long time to complete one instance due to its complexity.
- Subdirectories containing more utilities

## Subdirectories:

`data` - Rotation curve data used for fitting

`Inputs` - ??????? Should probably remove this, I think I made it but don't remember why

`NotFromTheLibrary` - Fitting attempts done without using custom NGC5533_functions library

`Static_Notebooks` - Notebooks which are not meant to be changed extensively. These exist for the sake of informing us which combinations of fitting parameters worked at a given point in time.

## Warnings

The halo.hdf5 file is already at maximum size; you will need to use `git checkout` instead of `git push` or you will run into errors.

In general, the .hdf5 files and git do not get along very well. For example, you cannot merge changes on an .hdf5 file, you will have to simply keep one or the other.