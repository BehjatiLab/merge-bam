from setuptools import setup, find_packages

setup(
    name='my_bam_tools',
    version='0.1',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    entry_points={
        'console_scripts': [
            'create-test-bed=create_position_bed:main',
            'merge-bam=merge_bam:main',
            'bed-from-txt=bed_from_mut_txt:main',
        ],
    },
    install_requires=[
        'pysam',  # Add other dependencies if necessary
    ],
    python_requires='>=3.11',
)
