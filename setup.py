from setuptools import setup

setup(
    name="dvartk",
    version="0.1.3c",
    description="variant comparison toolkit",
    url="http://github.com/soymintc/dvartk",
    author="Seongmin Choi",
    author_email="soymintc@gmail.com",
    license="MIT",
    packages=["dvartk"],
    install_requires=[
        "pandas",
        "numpy",
        "matplotlib",
        "seaborn",
        "pyfaidx",
    ],
    zip_safe=False,
)
