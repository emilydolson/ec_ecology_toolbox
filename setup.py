
# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

__version__ = "0.0.8"

# The main interface is through Pybind11Extension.
# * You can add cxx_std=11/14/17, and then build_ext can be removed.
# * You can set include_pybind11=false to add the include directory yourself,
#   say from a submodule.
#
# Note:
#   Sort input source files if you glob sources to ensure bit-for-bit
#   reproducible builds (https://github.com/pybind/python_example/pull/53)

ext_modules = [
    Pybind11Extension("ec_ecology_toolbox.selection_probabilities",
                      ["source/wrapper.cc"],
                      include_dirs=['third-party/Empirical/include/emp'],
                      cxx_std=20  # ,
                    #   define_macros=[('_DEBUG', 1), ('DEBUG', 1), ('NDEBUG', 0)],
                    #   extra_compile_args=["-g", "-O0"]
                      ),
]

setup(
    name="ec_ecology_toolbox",
    version=__version__,
    author="Emily Dolson",
    author_email="emilyldolson@gmail.com",
    url="https://github.com/emilydolson/ec_ecology_toolbox",
    description="Tools to analyze the ecology of evolutionary algorithms",
    long_description="Tools to analyze the ecology of evolutionary algorithms",
    ext_modules=ext_modules,
    packages=['ec_ecology_toolbox', 'ec_ecology_toolbox.community_assembly_graph'],

    extras_require={"test": "pytest"},
    # Currently, build_ext only provides an optional "highest supported C++
    # level" feature, but in the future it may provide more features.
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.7"
)
