# -*- coding: utf-8 -*-


from setuptools import setup


setup(
    name="pvutils",
    version="0.1",
    author="Matthias Mayr, Ivo Steinbrecher",
    author_email="matthias.mayr@unibw.de, ivo.steinbrecher@unibw.de",
    description="Utils for ParaView",
    packages=["blender", "pvutils"],
    install_requires=[
        "numpy==1.21.1",  # This is the version provided with ParaView 5.10.1
    ],
)
