#ifndef _INTERACTION_NETWORKS_H_
#define _INTERACTION_NETWORKS_H_

#include <functional>
#include <map>
#include <limits>
#include <unordered_set>

#include "base/vector.hpp"
#include "base/assert.hpp"
#include "datastructs/vector_utils.hpp"
#include "datastructs/map_utils.hpp"
#include "tools/string_utils.hpp"
#include "datastructs/set_utils.hpp"
#include "math/math.hpp"
#include "math/distances.hpp"
// #include "datastructs/Graph.hpp"
#include "Evolve/Resource.hpp"

#include <math.h>  

// from https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template<class T, class U>
typename std::enable_if<!std::numeric_limits<T>::is_integer || !std::numeric_limits<U>::is_integer, bool>::type
    almost_equal(T x, U y, int ulp)
{
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x-y) <= std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
    // unless the result is subnormal
           || std::abs(x-y) < std::numeric_limits<T>::min();
}

// integer version - just use the equals operator
template<class T, class U>
typename std::enable_if<std::numeric_limits<T>::is_integer && std::numeric_limits<U>::is_integer, bool>::type
    almost_equal(T x, U y, int ulp)
{
    return x == y;
}

// from https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template<class T, class U>
typename std::enable_if<!std::numeric_limits<T>::is_integer || !std::numeric_limits<U>::is_integer, bool>::type
    lt_almost_equal(T x, U y, int ulp)
{
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return x-y <= std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
    // unless the result is subnormal
           || x-y < std::numeric_limits<T>::min();
}

// integer version - just use the equals operator
template<class T, class U>
typename std::enable_if<std::numeric_limits<T>::is_integer && std::numeric_limits<U>::is_integer, bool>::type
    lt_almost_equal(T x, U y, int ulp)
{
    return x <= y;
}

// from https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template<class T, class U>
typename std::enable_if<!std::numeric_limits<T>::is_integer || !std::numeric_limits<U>::is_integer, bool>::type
    gt_almost_equal(T x, U y, int ulp)
{
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return y-x <= std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
    // unless the result is subnormal
           || y-x < std::numeric_limits<T>::min();
}

// integer version - just use the equals operator
template<class T, class U>
typename std::enable_if<std::numeric_limits<T>::is_integer && std::numeric_limits<U>::is_integer, bool>::type
    gt_almost_equal(T x, U y, int ulp)
{
    return x >= y;
}

/// Find the elements in @param pop with the highest value on
/// @param axis. PHEN_T should be a vector of some numeric type.
/// If epsilon lexicase is being used, specify a value of @param epsilon
/// to include everything within that range of the highest value
template <typename PHEN_T>
emp::vector<int> FindHighestIndices(emp::vector<PHEN_T> & pop, int axis, double epsilon = 0) {
    double best = std::numeric_limits<double>::lowest();
    emp::vector<int> winners;

    for (size_t i = 0; i < pop.size(); i++) {
        if (pop[i][axis] > best) {
            best = pop[i][axis];
            winners.resize(0);
            winners.push_back(i);
        } else if (almost_equal(pop[i][axis], best, 5)) {
            winners.push_back(i);
        }
    }

    if (epsilon) {
        winners.resize(0);
        for (size_t i = 0; i < pop.size(); i++) {
            if (pop[i][axis] + epsilon > best) {
                winners.push_back(i);
            }
        }
    }
    return winners;
}

template <typename PHEN_T>
bool IsElite(emp::vector<PHEN_T> & pop, int axis, PHEN_T individual, double epsilon = 0) {
    double best = individual[axis];

    for (size_t i = 0; i < pop.size(); i++) {
        if (pop[i][axis] > best + epsilon) {
            return false;
        }
    }
    return true;
}

template <typename PHEN_T>
emp::vector<int> FindWinningAxes(emp::vector<PHEN_T> & pop, emp::vector<int> & axes, PHEN_T individual, double epsilon = 0) {
    emp::vector<int> winning;
    for (int ax : axes) {
        if (IsElite(pop, ax, individual, epsilon)) {
            winning.push_back(ax);
        }
    }
    return winning;
}

template <typename PHEN_T>
emp::vector<PHEN_T> FindHighest(emp::vector<PHEN_T> & pop, int axis, double epsilon = 0) {
    // From https://stackoverflow.com/questions/9650679/c-select-a-subset-of-a-stdvector-based-predefined-element-indices

    emp::vector<int> highest = FindHighestIndices(pop, axis, epsilon);
    std::vector<PHEN_T> result(highest.size());

    std::transform(highest.begin(), highest.end(), result.begin(), [&pop](size_t pos) {return pop[pos];});
    return result;
}

/// Multiply all of the ints in @param nums together.
long int VectorProduct(const emp::vector<int> & nums) {
    long int result = 1;

    for (int num : nums) {
        result *= num;
    }

    return result;
}

template <typename PHEN_T>
void FilterImpossible(emp::vector<PHEN_T> & pop, emp::vector<int> & axes, double epsilon = 0) { 
    emp::vector<int> best;
    emp::vector<int> temp;
    for (int ax : axes) {
        temp = FindHighestIndices(pop, ax, epsilon);
        best.reserve(best.size() + distance(temp.begin(),temp.end()));
        best.insert(best.end(),temp.begin(),temp.end());       
    }

    best = emp::RemoveDuplicates(best);
    std::vector<PHEN_T> result(best.size());

    std::transform(best.begin(), best.end(), result.begin(), [&pop](size_t pos) {return pop[pos];});

    pop = result;
}

template <typename PHEN_T>
void FilterDominated(emp::vector<PHEN_T> & pop, emp::vector<int> & axes, double epsilon = 0) {
    emp::vector<int> to_remove;
    
    for (int i = 0; i < pop.size(); i++) {
        for (int j = 0; j < pop.size(); j++) {
            if (i == j) {
                continue;
            }
            bool losses = false;
            bool wins = false;
            for (int ax : axes) {
                if (pop[i][ax] + epsilon < pop[j][ax]) {
                    losses = true;
                } else if (pop[i][ax] - epsilon > pop[j][ax]) {
                    wins = true;
                    break;
                } 
            }
            if (losses && !wins) {
                to_remove.push_back(i);
                break;
            }
        }
    }
    // if (pop.size() < 10 && axes.size() < 6) {
    //     std::cout << emp::to_string(pop) << emp::to_string(axes) << emp::to_string(to_remove) << std::endl;
    // }
    for (emp::vector<int>::reverse_iterator i = to_remove.rbegin(); i != to_remove.rend(); ++i ) {
        pop.erase(pop.begin()+(*i));
    }
}


template <typename PHEN_T>
void PruneAxes(emp::vector<int> & axes, emp::vector<PHEN_T> & pop, double epsilon = 0) {
    
    emp::vector<int> to_remove;

    for (int ax : axes) {
        // std::cout << "Evaluating " << ax << std::endl;
        bool include = false;
        double best = pop[0][ax];
        double lowest = pop[0][ax];
        for (PHEN_T & org : pop) {
            if (org[ax] > best) {
                // std::cout << org[ax] << " > "<< best << std::endl;
                if (org[ax] - epsilon > lowest) {
                    // std::cout <<" Including  "<< ax << std::endl;
                    include = true;
                    break;
                }
                best = org[ax];
            } else if (org[ax] < lowest) {
                // std::cout << org[ax] << " < "<< lowest << std::endl;
                if (org[ax] + epsilon < best) {
                    // std::cout <<" Including  "<< ax << std::endl;
                    include = true;
                    break;
                }
                lowest = org[ax];
            }
        }

        if (!include) {
            to_remove.push_back(ax);
        }
    }

    for (int ax : to_remove) {
        axes.erase(std::remove(axes.begin(), axes.end(), ax), axes.end());
    }
}

template <typename PHEN_T>
void DeDuplicateAxes(emp::vector<int> & axes, emp::vector<PHEN_T> & pop, std::map<int, int> & dups) {
    emp::vector<int> to_remove;
    for (int ax : axes) {
        if (emp::Has(to_remove, ax)) {
            continue;
        }
        for (int ax2 : axes) {
            if (emp::Has(to_remove, ax2)) {
                continue;
            }
            if (ax == ax2) {
                continue;
            }
            bool same = true;
            for (PHEN_T & org : pop) {
                if (org[ax] != org[ax2]) {
                    same = false;
                    break;
                }
            }
            if (same) {
                to_remove.push_back(ax2);
                dups[ax] += dups[ax2];
            }
        }
    }

    for (int ax : to_remove) {
        axes.erase(std::remove(axes.begin(), axes.end(), ax), axes.end());
    }

}

template <typename PHEN_T>
void HandleTwoOrgs(std::map<PHEN_T, double> & fit_map, emp::vector<PHEN_T> & winners, emp::vector<int> axes, emp::vector<int> perm_levels, std::map<int,int> dups, double epsilon = 0, double multiplier=1.0) {
    double wins = 0;
    double losses = 0;
    for (int ax : axes) {
        if (winners[0][ax] > winners[1][ax] + epsilon) {
            wins += dups[ax];
        } else if (winners[1][ax] > winners[0][ax] + epsilon) {
            losses += dups[ax];
        }
    }
    if (wins > 0 || losses > 0) {
        fit_map[winners[0]] += (wins/(wins+losses)) * (multiplier/VectorProduct(perm_levels));//(double)emp::Factorial(axes.size() - 1);
        fit_map[winners[1]] += (losses/(wins+losses)) * (multiplier/VectorProduct(perm_levels));//(double)emp::Factorial(axes.size() - 1);
    } else {
        fit_map[winners[0]] += .5 * (multiplier/VectorProduct(perm_levels));//(double)emp::Factorial(axes.size() - 1);
        fit_map[winners[1]] += .5 * (multiplier/VectorProduct(perm_levels));//(double)emp::Factorial(axes.size() - 1);                
    }
}

template <typename PHEN_T>
double HandleTwoOrgsIndividual(emp::vector<PHEN_T> & winners, emp::vector<int> axes, emp::vector<int> perm_levels, double epsilon = 0) {
    double wins = 0;
    double losses = 0;

    for (int ax : axes) {
        if (winners[0][ax] > winners[1][ax] + epsilon) {
            wins++;
        } else if (winners[1][ax] > winners[0][ax] + epsilon) {
            losses++;
        }
    }
    if (wins > 0 || losses > 0) {
        return (wins/(wins+losses)) * (1.0/VectorProduct(perm_levels));//(double)emp::Factorial(axes.size() - 1);
    } else {
        return .5 * (1.0/VectorProduct(perm_levels));//(double)emp::Factorial(axes.size() - 1);
    }
}

template <typename PHEN_T>
double TraverseDecisionTreeIndividual(emp::vector<PHEN_T> & pop, emp::vector<int> axes, PHEN_T individual, emp::vector<int> perm_levels, double epsilon = 0) {
    // std::cout << "begining round of recursion " << axes.size() << emp::to_string(pop) << emp::to_string(axes) << std::endl;

    // std::cout << emp::to_string(pop) << emp::to_string(axes) << std::endl;
 
    emp_assert(pop.size() > 0, axes, perm_levels);

    // There's only one fitness criterion left, so it wins
    if (axes.size() == 1) {
        emp::vector<PHEN_T> winners = FindHighest(pop, axes[0], epsilon);
        if (emp::Has(winners, individual)) {
            return 1.0/((double)winners.size()*VectorProduct(perm_levels));
        }
        return 0;
    }

    // There's only one thing in the population, so it wins
    if (pop.size() == 1) {
        if (pop[0] == individual) {
            return 1.0/VectorProduct(perm_levels);
        }
        return 0;
    }

    perm_levels.push_back(axes.size());

    double plex = 0;
    for (int ax : FindWinningAxes(pop, axes, individual, epsilon)) {
        emp::vector<PHEN_T> winners = FindHighest(pop, ax, epsilon);
        emp::vector<int> next_axes = axes;
        emp_assert(emp::Has(winners, individual), pop, winners, individual, epsilon, ax);
        next_axes.erase(std::remove(next_axes.begin(), next_axes.end(), ax), next_axes.end());
        FilterImpossible(winners, next_axes, epsilon);

        // Check whether individual got filtered out by FilterImpossible
        // which can happen if ax was the last remaining criterion it wins on
        if (!emp::Has(winners, individual)) {
            continue;
        }

        if (winners.size() == 1) { // Not a tie
            plex += 1.0/VectorProduct(perm_levels);

        } else if (winners.size() == 2) { // optimization
            if (winners[0] != individual) {
                winners[1] = winners[0];
                winners[0] = individual;
            }
            plex += HandleTwoOrgsIndividual(winners, next_axes, perm_levels, epsilon);

        } else { // tie
            plex += TraverseDecisionTreeIndividual(winners, next_axes, individual, perm_levels, epsilon);
        }
    }
    return plex;
}

template <typename PHEN_T>
void TraverseDecisionTree(std::map<PHEN_T, double> & fit_map, emp::vector<PHEN_T> & pop, emp::vector<int> axes, emp::vector<int> perm_levels, std::map<int,int> dups, double epsilon = 0, double multiplier = 1.0) {
    // std::cout << "begining round of recursion " << axes.size() << emp::to_string(pop) << emp::to_string(axes) << std::endl;

    // std::cout << emp::to_string(pop) << emp::to_string(axes) << std::endl;
 
    emp_assert(pop.size() > 0, axes, perm_levels);
    emp_assert(axes.size() > 0, axes, perm_levels);

    // There's only one fitness criterion left, so it wins
    if (axes.size() == 1) {
        emp::vector<PHEN_T> winners = FindHighest(pop, axes[0], epsilon);
        // std::cout << "Winners: " << emp::to_string(winners) << std::endl;
        for (PHEN_T & org : winners) {
            fit_map[org]+=multiplier*dups[axes[0]]/((double)winners.size()*VectorProduct(perm_levels));
        }
        return;
    }

    // There's only one thing in the population, so it wins
    if (pop.size() == 1) {
        // TODO: Figure out if we actually need this
        // emp::vector<PHEN_T> winners = FindHighest(pop, axes[0], epsilon);
        // std::cout << "Winners: " << emp::to_string(winners) << std::endl;
        for (PHEN_T & org : pop) {
            fit_map[org]+=multiplier/VectorProduct(perm_levels);
        }
        return;
    }

    PruneAxes(axes, pop, epsilon);
    DeDuplicateAxes(axes, pop, dups);

    int real_axes = 0;
    for (int ax : axes){
        real_axes += dups[ax];
    }

    if (axes.size() == 0) {
        for (PHEN_T & org : pop) {
            fit_map[org]+=multiplier/((double)pop.size()*VectorProduct(perm_levels));
        }        
        return;
    }

    perm_levels.push_back(real_axes);

    // Check for only one axis again now that we've pruned the axes
    if (axes.size() == 1) {
        emp::vector<PHEN_T> winners = FindHighest(pop, axes[0], epsilon);
        // std::cout << "Winners: " << emp::to_string(winners) << std::endl;
        for (PHEN_T & org : winners) {
            fit_map[org]+=multiplier*real_axes*dups[axes[0]]/((double)winners.size()*VectorProduct(perm_levels));
        }
        return;
    }

    // std::cout << "post processing: " << axes.size() << emp::to_string(pop) << emp::to_string(axes) << std::endl;

    for (int ax : axes) {
        // if (perm_levels.size() == 1) {
        //     std::cout << "Axis: " << ax << " out of " << emp::to_string(axes) << std::endl;
        // }
        emp::vector<PHEN_T> winners = FindHighest(pop, ax, epsilon);
        emp::vector<int> next_axes = axes;
        next_axes.erase(std::remove(next_axes.begin(), next_axes.end(), ax), next_axes.end());
        if (epsilon) {
            FilterImpossible(winners, next_axes, epsilon);
        } else {
            FilterDominated(winners, next_axes, epsilon);
        }
        
        double new_multiplier = multiplier * dups[ax];

        if (winners.size() == 1) { // Not a tie
            // std::cout << "1 winner: " << emp::to_string(winners[0]) << " Controls " << (double)emp::Factorial(axes.size() - 1)<< std::endl;
            fit_map[winners[0]] += new_multiplier/VectorProduct(perm_levels);//(double)emp::Factorial(axes.size() - 1);
        } else if (winners.size() == 2) { // optimization
            HandleTwoOrgs(fit_map, winners, next_axes, perm_levels, dups, epsilon, new_multiplier);
        // } else if (winners.size() < 100 && next_axes.size() > 2) { // optimization
        //     HandleThreeOrgs(fit_map, winners, next_axes, perm_levels, epsilon);
        } else { // tie
            TraverseDecisionTree(fit_map, winners, next_axes, perm_levels, dups, epsilon, new_multiplier);
        }
    }
}

template <typename PHEN_T>
void TraverseDecisionTreeNaive(std::map<PHEN_T, double> & fit_map, emp::vector<PHEN_T> & pop, emp::vector<int> axes, emp::vector<int> perm_levels, double epsilon = 0) {
    // std::cout << "begining round of recursion " << axes.size() << emp::to_string(pop) << emp::to_string(axes) << std::endl;

    // std::cout << emp::to_string(pop) << emp::to_string(axes) << std::endl;
 
    emp_assert(pop.size() > 0, axes, perm_levels);

    if (axes.size() == 0) {
        for (PHEN_T & org : pop) {
            fit_map[org]+=1.0/(pop.size() * VectorProduct(perm_levels));
        }
        return;        
    }

    // There's only one thing in the population, so it wins
    if (pop.size() == 1) {
        emp::vector<PHEN_T> winners = FindHighest(pop, axes[0], epsilon);
        // std::cout << "Winners: " << emp::to_string(winners) << std::endl;
        for (PHEN_T & org : winners) {
            fit_map[org]+=1.0/VectorProduct(perm_levels);
        }
        return;
    }

    perm_levels.push_back(axes.size());

    for (int ax : axes) {
        emp::vector<PHEN_T> winners = FindHighest(pop, ax, epsilon);
        emp::vector<int> next_axes = axes;
        next_axes.erase(std::remove(next_axes.begin(), next_axes.end(), ax), next_axes.end());
        TraverseDecisionTreeNaive(fit_map, winners, next_axes, perm_levels, epsilon);
    }
}



template <typename PHEN_T>
emp::vector<double> LexicaseFitness(emp::vector<PHEN_T> & pop, double epsilon = 0) {

    if (pop.size() == 0) {
        return emp::vector<double>();
    }
    std::map<PHEN_T, double> fit_map;
    size_t n_funs = pop[0].size();

    for (PHEN_T & org : pop) {
        fit_map[org] = 0.0;
    }

    std::map<int,int> dups;
    for (int i = 0; i < n_funs; i++) {
        dups[i] = 1;
    }

    emp::vector<PHEN_T> de_dup_pop = emp::RemoveDuplicates(pop);
    TraverseDecisionTree(fit_map, de_dup_pop, emp::NRange(0, (int)n_funs), {}, dups, epsilon);

    for (PHEN_T & org : de_dup_pop) {
        fit_map[org] /= emp::Count(pop, org);
    }

    emp::vector<double> result;
    for (PHEN_T & org : pop) {
        result.push_back(fit_map[org]);
    }

    return result;
}

template <typename PHEN_T>
emp::vector<double> UnoptimizedLexicaseFitness(emp::vector<PHEN_T> & pop, double epsilon = 0) {

    emp_assert(pop.size() > 0);
    std::map<PHEN_T, double> fit_map;
    size_t n_funs = pop[0].size();

    for (PHEN_T & org : pop) {
        fit_map[org] = 0.0;
    }

    emp::vector<PHEN_T> de_dup_pop = emp::RemoveDuplicates(pop);
    TraverseDecisionTreeNaive(fit_map, de_dup_pop, emp::NRange(0, (int)n_funs), {}, epsilon);

    for (PHEN_T & org : de_dup_pop) {
        fit_map[org] /= emp::Count(pop, org);
    }

    emp::vector<double> result;
    for (PHEN_T & org : pop) {
        result.push_back(fit_map[org]);
    }

    return result;
}

template <typename PHEN_T>
double LexicaseFitnessIndividual(emp::vector<PHEN_T> & pop, int i, double epsilon = 0) {

    emp_assert(pop.size() > 0);
    size_t n_funs = pop[0].size();

    PHEN_T phen = pop[i];
    double n_dups = emp::Count(pop, phen);
    emp::vector<PHEN_T> de_dup_pop = emp::RemoveDuplicates(pop);
    return TraverseDecisionTreeIndividual(de_dup_pop, emp::NRange(0, (int)n_funs), phen, {}, epsilon)/n_dups;
}


void TournamentHelper(emp::vector<double> & fit_map, int t_size = 2){

    emp::vector<double> base_fit_map = fit_map;
    double pop_size = base_fit_map.size();

    for (size_t i = 0; i < fit_map.size(); i++) {
        double less = 0.0;
        double equal = 0.0; 

        for (auto & org2 : base_fit_map) {
            if (almost_equal(org2, base_fit_map[i], 5)) {
                equal++;
            } else if (org2 < base_fit_map[i]) {
                less++;
            }
        }

        long double p_less = less/pop_size;
        long double p_equal = equal/pop_size;
        long double p_self = 1/pop_size;

        fit_map[i] = (pow(p_equal + p_less, t_size) - pow(p_less, t_size)) * p_self/p_equal;
    }
}

template <typename PHEN_T>
emp::vector<double> RandomFitness(emp::vector<PHEN_T> & pop) {
    emp_assert(pop.size() > 0);
    emp::vector<double> fit_map;
    for (PHEN_T & org : pop) {
        fit_map.push_back(1.0/(double)pop.size());
    }
    // TODO: Double check
    TournamentHelper(fit_map, 1);
    return fit_map;
}

template <typename PHEN_T>
emp::vector<double> TournamentFitness(emp::vector<PHEN_T> & pop, int t_size = 2) {
    emp_assert(pop.size() > 0);
    emp::vector<double> fit_map;
    for (PHEN_T & org : pop) {
        fit_map.push_back(emp::Sum(org));
    }
    TournamentHelper(fit_map, t_size);
    return fit_map;
}

template <typename PHEN_T>
emp::vector<double> SharingFitness(emp::vector<PHEN_T> & pop, int t_size = 2, double alpha = 1, double sigma_share = 8.0) {
    // std::cout << "SHARING" << std::endl;
    emp::vector<double> fit_map;

    // std::cout << SigmaShare::Get(attrs) << std::endl;

    for (size_t i = 0; i < pop.size(); i++) {
        fit_map.push_back(1.0);
    }

    for (size_t i = 0; i < pop.size(); i++) {
        double niche_count = 0;

        for (PHEN_T & org2 : pop) {
            // Sharing function is euclidean distance
            // we could make this configurable
            double dist = emp::EuclideanDistance(pop[i], org2);
            if (dist < sigma_share) {
                niche_count += 1 - pow((dist/sigma_share), alpha);
            } 
        }

        // Don't worry, niche_count will always be at least one because each individual
        // increases its own niche count by 1
        fit_map[i] = emp::Sum(pop[i]) / niche_count;
    }
    TournamentHelper(fit_map, t_size);

    return fit_map;    
};


// TODO: the handling of resource use here needs to be re-thought
// template <typename PHEN_T>
// emp::vector<double> EcoEAFitness(emp::vector<PHEN_T> & pop, all_attrs attrs=DEFAULT) {
//     // std::cout << "ECO-EA" << std::endl;
//     emp_assert(pop.size() > 0);

//     emp::vector<double> fit_map;

//     size_t n_funs = pop[0].size();

//     for (size_t i = 0; i < pop.size(); i++) {
//         fit_map.push_back(1.0);
//     }

//     for (int axis : emp::NRange(0, (int)n_funs)) {
//         double res = ResourceInflow::Get(attrs)/(double)pop.size(); // Average resource available per individual
//         double count = 0;

//         for (PHEN_T & org : pop) {
//             if (org[axis] >= NicheWidth::Get(attrs)) {
//                 count+= std::min(Cf::Get(attrs)*res*pow(org[axis]/MaxScore::Get(attrs),2.0) - Cost::Get(attrs), MaxBonus::Get(attrs));
//             }
//         }

//         count /= (double) pop.size(); // average resource use per pop member
//         res -= count; // how much of average available resources are used on average
//         res = std::max(res, 0.0); // Can't have a negative amount of a resource

//         for (size_t i = 0; i < pop.size(); i++) {
//             if (pop[i][axis] >= NicheWidth::Get(attrs)) {
//                 fit_map[i] *= pow(2,std::min(Cf::Get(attrs)*res*pow(pop[i][axis]/MaxScore::Get(attrs),2.0) - Cost::Get(attrs), MaxBonus::Get(attrs)));
//             }
//         }
//     }

//     for (size_t i = 0 ; i < pop.size(); i++) {
//         std::cout << i << " " << fit_map[i] << std::endl;
//     }

//     TournamentHelper(fit_map, TournamentSize::Get(attrs));

//     return fit_map;
// };

template <typename PHEN_T>
std::function<emp::vector<double>(emp::vector<PHEN_T>&)> do_lexicase = [](emp::vector<PHEN_T> & pop){return LexicaseFitness(pop);};

template <typename PHEN_T>
std::function<emp::vector<double>(emp::vector<PHEN_T>&)> do_tournament = [](emp::vector<PHEN_T> & pop){return TournamentFitness(pop);};

// template <typename PHEN_T>
// std::function<emp::vector<double>(emp::vector<PHEN_T>&, all_attrs)> do_eco_ea = [](emp::vector<PHEN_T> & pop, all_attrs attrs=DEFAULT){return EcoEAFitness(pop,attrs);};

template <typename PHEN_T>
std::function<emp::vector<double>(emp::vector<PHEN_T>&)> do_sharing = [](emp::vector<PHEN_T> & pop){return SharingFitness(pop);};

template <typename PHEN_T>
std::function<emp::vector<double>(emp::vector<PHEN_T>&)> do_random = [](emp::vector<PHEN_T> & pop){return RandomFitness(pop);};


// template <typename PHEN_T>
// emp::WeightedGraph CalcCompetition(emp::vector<PHEN_T> pop, 
//                         std::function<emp::vector<double>(emp::vector<PHEN_T>&)> fit_fun) {

//     emp::WeightedGraph effects(pop.size());


//     emp::vector<double> fitnesses = fit_fun(pop);

//     for (size_t i = 0; i < pop.size(); i++) {
//         effects.SetLabel(i, emp::to_string(pop[i]));
//         // std::cout << effects.GetLabel(i) << std::endl;

//         emp::vector<PHEN_T> curr = pop;
//         for (size_t ax = 0; ax < curr[i].size(); ax++) {
//             curr[i][ax] = 0; // Replace org with null org so pop size stays same
//         }

//         emp::vector<double> new_fits = fit_fun(curr);
//         for (size_t j = 0; j < pop.size(); j++ ) {
//             if (i == j) {continue;}

//             // In terms of floating-point imprecision issues, it's much better
//             // to check for equality than doing the subtraction and checking
//             // for the result to be 0
//             if (!almost_equal(fitnesses[j], new_fits[j], 5)) {
//                 effects.AddEdge(i, j, fitnesses[j] - new_fits[j]);
//             }
//             // std::cout << effect << std::endl;
//         }
//     }

//     return effects;
// }


#endif
