#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "../Catch2/single_include/catch2/catch.hpp"
// #include "../Catch2/include/internal/benchmark/catch_benchmark.hpp"

#include "../source/nsga_helpers.hpp"
#include "datastructs/vector_utils.hpp"

TEST_CASE("Test NSGA") {
    emp::vector<emp::vector<double>> pop({{1,0,0}, 
                                          {0,1,0},
                                          {0,0,1},
                                          {0,0,0}});

    emp::vector<double> fit_map(4);
    emp::vector<emp::vector<double>> result = NSGAIIHelper(fit_map, pop, 3);
    CHECK(emp::Has(result, emp::vector<double>({1,0,0})));
    CHECK(emp::Has(result, emp::vector<double>({0,1,0})));
    CHECK(emp::Has(result, emp::vector<double>({0,0,1})));
    CHECK(!emp::Has(result, emp::vector<double>({0,0,0})));
    CHECK(fit_map[0] == 1);
    CHECK(fit_map[1] == 1);
    CHECK(fit_map[2] == 1);
    CHECK(fit_map[3] == 0);

    fit_map.clear();
    fit_map.resize(pop.size());
    result = NSGAIIHelper(fit_map, pop, 4);
    CHECK(emp::Has(result, emp::vector<double>({1,0,0})));
    CHECK(emp::Has(result, emp::vector<double>({0,1,0})));
    CHECK(emp::Has(result, emp::vector<double>({0,0,1})));
    CHECK(emp::Has(result, emp::vector<double>({0,0,0})));
    CHECK(fit_map[0] == 1);
    CHECK(fit_map[1] == 1);
    CHECK(fit_map[2] == 1);
    CHECK(fit_map[3] == 1);

    pop = emp::vector<emp::vector<double>>({{1,0,0}, 
                                            {0,1,0},
                                            {0,0,1},
                                            {2,0,0}, 
                                            {0,2,0},
                                            {0,0,2},
                                            {0,0,0}});

    fit_map.clear();
    fit_map.resize(pop.size());
    result = NSGAIIHelper(fit_map, pop, 3);
    CHECK(emp::Has(result, emp::vector<double>({2,0,0})));
    CHECK(emp::Has(result, emp::vector<double>({0,2,0})));
    CHECK(emp::Has(result, emp::vector<double>({0,0,2})));
    CHECK(!emp::Has(result, emp::vector<double>({1,0,0})));
    CHECK(!emp::Has(result, emp::vector<double>({0,1,0})));
    CHECK(!emp::Has(result, emp::vector<double>({0,0,1})));
    CHECK(!emp::Has(result, emp::vector<double>({0,0,0})));
    CHECK(fit_map[0] == 0);
    CHECK(fit_map[1] == 0);
    CHECK(fit_map[2] == 0);
    CHECK(fit_map[3] == 1);
    CHECK(fit_map[4] == 1);
    CHECK(fit_map[5] == 1);
    CHECK(fit_map[6] == 0);

    fit_map.clear();
    fit_map.resize(pop.size());
    result = NSGAIIHelper(fit_map, pop, 6);
    CHECK(emp::Has(result, emp::vector<double>({2,0,0})));
    CHECK(emp::Has(result, emp::vector<double>({0,2,0})));
    CHECK(emp::Has(result, emp::vector<double>({0,0,2})));
    CHECK(emp::Has(result, emp::vector<double>({1,0,0})));
    CHECK(emp::Has(result, emp::vector<double>({0,1,0})));
    CHECK(emp::Has(result, emp::vector<double>({0,0,1})));
    CHECK(!emp::Has(result, emp::vector<double>({0,0,0})));
    CHECK(fit_map[0] == 1);
    CHECK(fit_map[1] == 1);
    CHECK(fit_map[2] == 1);
    CHECK(fit_map[3] == 1);
    CHECK(fit_map[4] == 1);
    CHECK(fit_map[5] == 1);
    CHECK(fit_map[6] == 0);

    pop = emp::vector<emp::vector<double>>({{1,1,2,0,0}, 
                                            {0,0,1,1,0},
                                            {0,0,0,2,1}});
    fit_map.clear();
    fit_map.resize(pop.size());
    result = NSGAIIHelper(fit_map, pop, 2);
    CHECK(emp::Has(result, emp::vector<double>({1,1,2,0,0})));
    CHECK(emp::Has(result, emp::vector<double>({0,0,0,2,1})));
    CHECK(!emp::Has(result, emp::vector<double>({0,0,1,1,0})));
    CHECK(fit_map[0] == 1);
    CHECK(fit_map[1] == 0);
    CHECK(fit_map[2] == 1);


}