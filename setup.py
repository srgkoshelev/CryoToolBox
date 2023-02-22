from setuptools import setup, find_packages

setup(
    name='heat_transfer',
    version='0.1.0',
    description='A Python package for heat transfer and fluid dynamics calculations',
    author='Sergey Koshelev',
    author_email='srg.koshelev@pm.me',
    url='https://github.com/srgkoshelev/heat_transfer',
    packages=find_packages('heat_transfer'),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
    ],
    install_requires=[
        'CoolProp',
        'Pint',
        'Serialize',
        'scipy',
        'pyyaml',
        'xlsxwriter',
    ],
)
