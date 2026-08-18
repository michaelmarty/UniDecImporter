"""Minimal HDF5 helpers used by :mod:`MassSpecImporter.FileParser`."""


def replace_dataset(group, name, data):
    """Replace *name* in an h5py group and return the new dataset."""
    if name in group:
        del group[name]
    return group.create_dataset(name, data=data)

