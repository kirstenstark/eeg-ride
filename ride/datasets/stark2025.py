import pooch

LOCAL_CACHE = 'eeg-ride'
BASE_URL = 'https://files.de-1.osf.io/v1/resources/289ej/providers/osfstorage'


def get_stark2025(path=None):
    """Get sample data from the Stark (2025) dataset.

    Data that are not yet available locally will be downloaded from the OSF.
    See :footcite:`???` for details on the dataset.

    Parameters
    ----------
    path : str or Path, optional
        Local directory path to download the data to. By default, uses the
        user's local cache directory. An alternative way to specify the
        download path is to set the environment variable ``RIDE_DATA_DIR``.

    Returns
    -------
    list
        A list with the file paths of the downloaded log file and the EEG
        header file.

    References
    ----------
    .. footbibliography::
    """

    if path is None:
        path = pooch.os_cache(LOCAL_CACHE)
        env = 'RIDE_DATA_DIR'
    else:
        env = None

    fetcher = pooch.create(path=path, base_url=BASE_URL, env=env)
    local_dir = fetcher.abspath

    local_paths = ['2_memory.txt',
                   'VP0302.vhdr',
                   'VP0302.vmrk',
                   'VP0302.eeg']
    hashes = ['md5:ef5f33e85f99d5344e208f67b2a6122e',
              'md5:38b92d3f68ebc045fe3c0ab6300720ac',
              'md5:e3bb2be657f1f836553fbb9ffe94be57',
              'md5:07836bf4854a1585ced36c07b3e1c1a5']
    urls = ['https://files.de-1.osf.io/v1/resources/289ej/providers/osfstorage/6797c7014e87d78691df2fd9',
            'https://files.de-1.osf.io/v1/resources/289ej/providers/osfstorage/6797c1006bf0c8d37bdf2d49',
            'https://files.de-1.osf.io/v1/resources/289ej/providers/osfstorage/6797c103794a3a633aa48b69',
            'https://files.de-1.osf.io/v1/resources/289ej/providers/osfstorage/6797c19d8099c9814c7c2316']
    local_files = []

    for local_path, hash, url in zip(local_paths, hashes, urls):

        local_file = local_dir.joinpath(local_path)

        if not local_file.exists():
            fetcher.registry[local_path] = hash
            fetcher.urls[local_path] = url
            _ = fetcher.fetch(local_path)

        local_files.append(local_file)

    return local_files[0:2]
