import pathlib

from abaqus_python.abaqus_interface import ABQInterface
from input_file_reader.input_file_reader import InputFileReader


abq = ABQInterface('abq2018')


def main():
    model_directory = pathlib.Path.home() / "python_projects" / "python_fatigue" / "planetary_gear" / "input_files"
    coarse_input_file = model_directory / 'coarse_tooth.inp'
    coarse_tooth = InputFileReader()
    coarse_tooth.read_input_file()