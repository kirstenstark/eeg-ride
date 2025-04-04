from setuptools import find_packages, setup

if __name__ == "__main__":

    # Paste README as long description
    with open('README.md', 'r', encoding='utf-8') as fh:
        long_description = fh.read()

    # Actual setup
    setup(
        name='eeg-ride',
        author='Kirsten Stark, Alexander Enge',
        author_email='kirsten.stark@hu-berlin.de, alexander.enge@hu-berlin.de',
        description='Removing speech artifacts from EEG data using Residue '
                    'Iteration Decomposition (RIDE)',
        long_description=long_description,
        long_description_content_type='text/markdown',
        url='https://github.com/kirstenstark/eeg-ride',
        project_urls={
            'Issue trackers': 'https://github.com/kirstenstark/eeg-ride/issues'
        },
        classifiers=[
            'Programming Language :: Python',
            'Programming Language :: Python :: 3',
            'License :: OSI Approved :: MIT License',
            'Operating System :: OS Independent',
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering'
        ],
        packages=find_packages(),
        install_requires=[
            'matplotlib',
            'mne',
            'numpy',
            'pandas',
            'scipy'
        ],
        python_requires='>=3.8',
    )
