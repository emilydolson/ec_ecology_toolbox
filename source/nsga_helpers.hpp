#ifndef PAGMO_NSGA_HELPERS_HPP
#define PAGMO_NSGA_HELPERS_HPP

/*From PAGMO - either rewrite or switch to GPL*/

#include "base/vector.hpp"

template <typename T, bool After = true>
bool less_than_f(T a, T b)
{
    static_assert(std::is_floating_point<T>::value, "less_than_f can be used only with floating-point types.");
    if (!std::isnan(a)) {
        if (!std::isnan(b))
            return a < b; // a < b
        else
            return After; // a < nan
    } else {
        if (!std::isnan(b))
            return !After; // nan < b
        else
            return false; // nan < nan
    }
}

template <typename T, bool After = true>
inline bool greater_than_f(T a, T b)
{
    static_assert(std::is_floating_point<T>::value, "greater_than_f can be used only with floating-point types.");
    if (!std::isnan(a)) {
        if (!std::isnan(b))
            return a > b; // a > b
        else
            return !After; // a > nan
    } else {
        if (!std::isnan(b))
            return After; // nan > b
        else
            return false; // nan > nan
    }
}

template <typename PHEN_T>
bool pareto_dominance(const PHEN_T &obj1, const PHEN_T &obj2)
{
    emp_assert(obj1.size() == obj2.size(), obj1, obj2);
    bool found_strictly_dominating_dimension = false;
    for (size_t i = 0; i < obj1.size(); ++i) {
        if (less_than_f(obj1[i], obj2[i])) {
            return false;
        } else if (greater_than_f(obj1[i], obj2[i])) {
            found_strictly_dominating_dimension = true;
        }
    }
    return found_strictly_dominating_dimension;
}

template <typename PHEN_T>
emp::vector<PHEN_T> crowding_distance(const emp::vector<PHEN_T> & non_dom_front)
{
    size_t N = non_dom_front.size();
    // We make sure to have two points at least
    size_t M = non_dom_front[0].size();
    // We make sure the first point of the input non dominated front contains at least two objectives

    emp::vector<size_t> indexes(N);
    std::iota(indexes.begin(), indexes.end(), 0);
    emp::vector<double> dists(N, 0.);
    for (size_t i = 0; i < M; ++i) {
        std::sort(indexes.begin(), indexes.end(), [i, &non_dom_front](size_t idx1, size_t idx2) {
            return greater_than_f(non_dom_front[idx1][i], non_dom_front[idx2][i]);
        });
        dists[indexes[0]] = std::numeric_limits<double>::infinity();
        dists[indexes[N - 1]] = std::numeric_limits<double>::infinity();
        double df = non_dom_front[indexes[N - 1]][i] - non_dom_front[indexes[0]][i];
        for (size_t j = 1; j < N - 1; ++j) {
            dists[indexes[j]] += (non_dom_front[indexes[j + 1]][i] - non_dom_front[indexes[j - 1]][i]) / df;
        }
    }

    std::sort(indexes.begin(), indexes.end(), [&dists](size_t idx1, size_t idx2) {
        return less_than_f(dists[idx1], dists[idx2]);
    }); // Ascending order

    emp::vector<PHEN_T> result;
    for (size_t i : indexes) {
        result.push_back(non_dom_front[i]);
    }

    return result;
}

template <typename PHEN_T>
emp::vector<PHEN_T> NSGAIIHelper(emp::vector<double> & fit_map, const emp::vector<PHEN_T> & pop, size_t N)
{
    emp::vector<PHEN_T> result;

    if (N == 0) {
        return result;
    } else if (pop.size() == 1 && N == 1) {
        result.push_back(pop[0]);
        return result;
    } else if (pop.size() == 2 && N==2) {
        if (pareto_dominance(pop[0], pop[1])) {
            result.push_back(pop[0]);
            return result;
        } else if (pareto_dominance(pop[1], pop[0])) {
            result.push_back(pop[1]);
            return result;
        } else {
            result.push_back(pop[0]);
            result.push_back(pop[1]);
            return result;
        }
    } 

    // Initialize the return values
    emp::vector<size_t> current_front;
    emp::vector<emp::vector<size_t>> dom_list(pop.size());
    emp::vector<size_t> dom_count(pop.size());

    // Start the fast non dominated sort algorithm
    for (size_t i = 0; i < pop.size(); ++i) {
        dom_list[i].clear();
        dom_count[i] = 0;
        for (size_t j = 0; j < i; ++j) {
            if (pareto_dominance(pop[i], pop[j])) {
                dom_list[i].push_back(j);
                ++dom_count[j];
            } else if (pareto_dominance(pop[j], pop[i])) {
                dom_list[j].push_back(i);
                ++dom_count[i];
            }
        }
    }

    for (size_t i = 0; i < pop.size(); ++i) {
        if (dom_count[i] == 0) {
            current_front.push_back(i);
        }
    }

    // If 0th front fits, add it to result
    if (current_front.size() <= N) {
        for (size_t i : current_front) {
            result.push_back(pop[i]);
        }
    } else {
        emp::vector<PHEN_T> front_phens;
        for (size_t i : current_front) {
            front_phens.push_back(pop[i]);
        }
        emp::vector<PHEN_T> sorted = crowding_distance(front_phens);
        while (result.size() < N) {
            result.push_back(sorted.back());
            sorted.pop_back();
        }
        return result;
    }

    while (current_front.size() != 0 && result.size() < N){
        std::vector<size_t> next_front;
        for (size_t p = 0; p < current_front.size(); ++p) {
            for (size_t q = 0; q < dom_list[current_front[p]].size(); ++q) {
                --dom_count[dom_list[current_front[p]][q]];
                if (dom_count[dom_list[current_front[p]][q]] == 0) {
                    next_front.push_back(dom_list[current_front[p]][q]);
                }
            }
        }

        current_front = next_front;
        if (result.size() + current_front.size() <= N) {
            for (size_t i : current_front) {
                result.push_back(pop[i]);
            }
        } else {
            emp::vector<PHEN_T> front_phens;
            for (size_t i : current_front) {
                front_phens.push_back(pop[i]);
            }
            emp::vector<PHEN_T> sorted = crowding_distance(front_phens);
            while (result.size() < N) {
                result.push_back(sorted.back());
                sorted.pop_back();
            }
        } 

    }
    return result;
}



#endif