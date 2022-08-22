from setuptools import setup, find_packages


setup(
    name="dvartk",
    version="0.1.8",
    license="MIT",
    author="Seongmin Choi",
    author_email="soymintc@gmail.com",
    packages=find_packages("src"),
    package_dir={"": "src"},
    url="https://github.com/soymintc/dvartk",
    keywords=["SNV", "SV", "comparison", "toolkit"],
    install_requires=[
        "pandas",
        "wgs-analysis",
        "numpy",
        "matplotlib",
        "matplotlib_venn",
        "seaborn",
    ],
)
