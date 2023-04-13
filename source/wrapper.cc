
#include "interaction_networks.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

PYBIND11_MODULE(ec_ecology_toolbox, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("LexicaseFitness", &LexicaseFitness<emp::vector<double>>, 
            R"mydelimiter(
            Return a vector containing the probability that each member of the population will be selected by lexicase selection or epsilon lexicase selection.
            
            For example, LexicaseFitness([[1, 2, 2], [2, 1, 2], [0, 0, 0]]) will return [.5, .5, 0], because the first two score vectors each have a 50% chance
            of being chosen and the third has no chance of being chosen.

            Note: calculating these probabilities is an NP-Hard problem (Dolson, 2023). This function is optimized, but if you try to use it with too large of input
            it might take a very a long time.

            Parameters
            ----------
            pop: A list of lists of floats indicating the scores of a population on each test case/fitness criterion.
            epsilon (optional): A floating point number indicating the epsilon value to use (if you want epsilon-lexicase selection probabilities; default value is 0, which is equivalent to standard lexicase selection).  
            )mydelimiter", 
          py::arg("pop"), py::arg("epsilon") = 0.0);
    // m.def("LexicaseFitness", [](emp::vector<emp::vector<double>> pop, double epsilon = 0){return LexicaseFitness(pop, epsilon);}, "Calculate probabilities of selection under lexicase selection");
}