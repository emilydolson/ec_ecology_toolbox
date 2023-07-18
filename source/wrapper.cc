
#include "interaction_networks.hpp"
#include "Evolve/NK.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

PYBIND11_MODULE(ec_ecology_toolbox, m) {
    m.doc() = "Tools for measuring the ecology of various evolutionary algorithms"; 

    m.def("LexicaseFitness", &LexicaseFitness<emp::vector<double>>, 
            R"mydelimiter(
            Return a vector containing the probability that each member of the population will be selected by lexicase selection or epsilon lexicase selection.
            
            For example, LexicaseFitness([[1, 2, 2], [2, 1, 2], [0, 0, 0]]) will return [.5, .5, 0], because the first two score vectors each have a 50% chance
            of being chosen and the third has no chance of being chosen.

            Note: calculating these probabilities is an NP-Hard problem (Dolson, 2023). This function is optimized, but if you try to use it with too large of input
            it might take a very a long time.

            Parameters
            ----------
            pop: list of lists of floats 
              The scores of each member of the population population on each test case/fitness criterion.
            epsilon: float
              (optional) The epsilon value to use (if you want epsilon-lexicase selection probabilities; default value is 0, which is equivalent to standard lexicase selection).  

            Returns
            -------
            List of floats
              The probabilities of each individual in pop being selected by lexicase selection.            
            )mydelimiter",
          py::arg("pop"), py::arg("epsilon") = 0.0);//, py::arg("binary") = false);
    
    m.def("UnoptimizedLexicaseFitness", &UnoptimizedLexicaseFitness<emp::vector<double>>, 
            R"mydelimiter(
            Return a vector containing the probability that each member of the population will be selected by lexicase selection or epsilon lexicase selection.
            
            For example, LexicaseFitness([[1, 2, 2], [2, 1, 2], [0, 0, 0]]) will return [.5, .5, 0], because the first two score vectors each have a 50% chance
            of being chosen and the third has no chance of being chosen.

            Note: calculating these probabilities is an NP-Hard problem (Dolson, 2023). This function is optimized, but if you try to use it with too large of input
            it might take a very a long time.

            Parameters
            ----------
            pop: list of lists of floats 
              The scores of each member of the population population on each test case/fitness criterion.
            epsilon: float
              (optional) The epsilon value to use (if you want epsilon-lexicase selection probabilities; default value is 0, which is equivalent to standard lexicase selection).  

            Returns
            -------
            List of floats
              The probabilities of each individual in pop being selected by lexicase selection.            
            )mydelimiter",
          py::arg("pop"), py::arg("epsilon") = 0.0);//, py::arg("binary") = false);

    m.def("LexicaseFitnessIndividual", &LexicaseFitnessIndividual<emp::vector<double>>, 
            R"mydelimiter(
            Returns the probability that a single individual is selected by lexicase selection.

            Note: calculating these probabilities is an NP-Hard problem (Dolson, 2023). This function is optimized, but if you try to use it with too large of input
            it might take a very a long time. This version is faster than LexicaseFitness if you just need the probability for a single individual, but is still worst-case O(N!)

            Parameters
            ----------
            pop: list of lists of floats 
              The scores of a population on each test case/fitness criterion.
            i: int
              The index of the individual in pop to calculate selection probability for
            epsilon: float
              (optional) The epsilon value to use (if you want epsilon-lexicase selection probabilities; default value is 0, which is equivalent to standard lexicase selection).  
            
            Returns
            -------
            float
              The probability of individual i being selected by lexicase selection.
            )mydelimiter", 
          py::arg("pop"), py::arg("i"), py::arg("epsilon") = 0.0);

        // m.def("LexicaseFitnessIndividualBinary", &SolveBinary<emp::vector<double>>, 
        //     R"mydelimiter(
        //     Returns the probability that a single individual is selected by lexicase selection in the case where all scores are either 0 or 1.

        //     Parameters
        //     ----------
        //     pop: list of lists of floats 
        //       The scores of a population on each test case/fitness criterion.
        //     i: int
        //       The index of the individual in pop to calculate selection probability for
            
        //     Returns
        //     -------
        //     float
        //       The probability of individual i being selected by lexicase selection.
        //     )mydelimiter", 
        //   py::arg("pop"), py::arg("i"));

    m.def("SharingFitness", &SharingFitness<emp::vector<double>>, 
            R"mydelimiter(
            Return a vector containing the probability that each member of the population will be selected under tournament selection with fitness sharing.

            The numbers in the pop parameter are assumed to be scores on a set of test cases/fitness criteria/tasks.
            Similarity will be calculated as the euclidean distance between these scores.
            Overall "Fitness" will be calculated as the sum of these scores, divided by the fitness sharing niche count.

            Parameters
            ----------
            pop: list of lists of floats 
              The scores of each member of the population population on each test case/fitness criterion.
            t_size: int
              Tournament size; the number of individuals that will be randomly selected to compete against each other in each selection event.
            alpha: float
              The alpha parameter of the fitness sharing function (controls shape of the sharing function)
            sigma_share: float
              The sharing threshold (i.e. how similar do individuals need to be to share fitness)

            Returns
            -------
            List of floats
              The probabilities of each individual in pop being selected.            
            )mydelimiter",
          py::arg("pop"), py::arg("t_size") = 2, py::arg("alpha") = 1, py::arg("sigma_share") = 8.0);

    m.def("TournamentFitness", &TournamentFitness<emp::vector<double>>, 
            R"mydelimiter(
            Return a vector containing the probability that each member of the population will be selected under tournament selection.

            The numbers in the pop parameter are assumed to be scores on a set of test cases/fitness criteria/tasks.
            Overall "Fitness" will be calculated as the sum of these scores.

            Parameters
            ----------
            pop: list of lists of floats 
              The scores of each member of the population population on each test case/fitness criterion.
            t_size: int
              Tournament size; the number of individuals that will be randomly selected to compete against each other in each selection event.

            Returns
            -------
            List of floats
              The probabilities of each individual in pop being selected.            
            )mydelimiter",
          py::arg("pop"), py::arg("t_size") = 2);

    py::class_<emp::Random>(m, "Random")
      .def(py::init<int>());

    py::class_<emp::BitVector>(m, "BitVector")
      .def(py::init<std::string>());

    py::class_<emp::NKLandscape>(m, "NKLandscape")
      .def(py::init<size_t, size_t, emp::Random&>())
      .def("GetFitness", static_cast<double (emp::NKLandscape::*)(size_t, size_t) const>(&emp::NKLandscape::GetFitness))
      .def("GetFitnesses", static_cast<emp::vector<double> (emp::NKLandscape::*)(emp::BitVector) const>(&emp::NKLandscape::GetFitnesses))
      ;
}