def get_data_filename(relative_path): #TODO put in utils
    """Get the full path to one of the reference files in testsystems.
    In the source distribution, these files are in ``mdfptools/data/``,
    but on installation, they're moved to somewhere in the user's python site-packages directory.

    Parameters
    ----------
    name : str
        Name of the file to load (with respect to the repex folder).

    Returns
    ---------
    fn : str
        filename
    """

    import os
    from pkg_resources import resource_filename
    fn = resource_filename('cpeptools', os.path.join('data', relative_path))

    if not os.path.exists(fn):
        raise ValueError("Sorry! %s does not exist. If you just added it, you'll have to re-install" % fn)

    return fn
