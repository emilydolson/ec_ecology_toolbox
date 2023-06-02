#ifndef _INTERACTION_NETWORKS_H_
#define _INTERACTION_NETWORKS_H_

#include <functional>
#include <map>
#include <limits>
#include <unordered_set>
#include <queue>

#include "base/vector.hpp"
#include "base/assert.hpp"
#include "datastructs/hash_utils.hpp"
#include "datastructs/vector_utils.hpp"
#include "datastructs/map_utils.hpp"
#include "tools/string_utils.hpp"
#include "datastructs/set_utils.hpp"
#include "math/math.hpp"
#include "math/distances.hpp"
#include "datastructs/Graph.hpp"
#include "Evolve/Resource.hpp"


#include <math.h>  

struct Axis {
    int dups = 1;
    int bonus_dups = 0;
    emp::vector<int> sub_axes;
    int orig_id;
    int curr_id;
    Axis(int id) {
        orig_id = id;
        curr_id = id;
    }
    bool operator==(const Axis & other) const {
        return curr_id == other.curr_id;
    }
    bool operator <(const Axis & other) const {
        return curr_id < other.curr_id;
    }
    friend std::ostream& operator<<(std::ostream& os, const Axis& a) {
        os << a.curr_id;
        return os;
    }
};

template <typename PHEN_T>
struct Individual {
    using value_type = typename PHEN_T::value_type;
    // PHEN_T orig_phen;
    PHEN_T curr_phen;
    int id;
    int dups = 1;

    Individual(PHEN_T p, int _id) {
        // orig_phen = p;
        curr_phen = p;
        id = _id;
    }
    Individual() {
        id = -1;
        curr_phen = {};
    }
    typename PHEN_T::value_type  operator[](int idx) const {return curr_phen[idx];}
    typename PHEN_T::value_type& operator[](int idx) {return curr_phen[idx];}
    bool operator==(const Individual & other) const {
        return curr_phen == other.curr_phen;
    }
    friend std::ostream& operator<<(std::ostream& os, const Individual& i) {
        os << i.curr_phen;
        return os;
    }
};

template<typename PHEN_T>
struct std::hash<Individual<PHEN_T>>
{
    std::size_t operator()(Individual<PHEN_T> const& s) const noexcept
    {
        return emp::ContainerHash<PHEN_T>{}(s.curr_phen);
    }
};



  /// Return a new vector containing the same elements as @param v
  /// with any duplicate elements removed.
  /// Not guaranteed to preserve order
  template <typename T>
  emp::vector<T> RemoveDuplicatesIndexes(const emp::vector<T> & v, emp::vector<Axis> & axes) {
    emp::vector<T> temp_vec = v;
    for (int i = 0; i < v.size(); i++) {
        temp_vec[i].curr_phen.resize(axes.size());
        for (int j = 0; j < axes.size(); j++) {
            temp_vec[i][j] = v[i][axes[j].curr_id];
        }
    }
    for (int i = 0; i < axes.size(); i++) {
        axes[i].curr_id = i;
    }
    std::unordered_set<T> temp_set(temp_vec.begin(), temp_vec.end());
    emp::vector<T> new_vec(temp_set.begin(), temp_set.end());
    std::cout << temp_set.size() << " " << temp_vec.size() << " " << temp_vec << " " << new_vec << std::endl;
    // emp_assert(new_vec.size() > 0, new_vec);
    return new_vec;
  }


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
            if (gt_almost_equal(pop[i][axis] + epsilon, best, 5)) {
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
        if (!lt_almost_equal(pop[i][axis], best + epsilon, 5)) {
            return false;
        }
    }
    return true;
}

template <typename PHEN_T>
emp::vector<Axis> FindWinningAxes(emp::vector<PHEN_T> & pop, emp::vector<Axis> & axes, PHEN_T individual, double epsilon = 0) {
    emp::vector<Axis> winning;
    for (Axis & ax : axes) {
        if (IsElite(pop, ax.curr_id, individual, epsilon)) {
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
void FilterImpossible(emp::vector<PHEN_T> & pop, emp::vector<Axis> & axes, double epsilon = 0) { 
    emp::vector<int> best;
    emp::vector<int> temp;
    for (Axis & ax : axes) {
        temp = FindHighestIndices(pop, ax.curr_id, epsilon);
        best.reserve(best.size() + distance(temp.begin(),temp.end()));
        best.insert(best.end(),temp.begin(),temp.end());       
    }

    best = emp::RemoveDuplicates(best);
    std::vector<PHEN_T> result(best.size());

    std::transform(best.begin(), best.end(), result.begin(), [&pop](size_t pos) {return pop[pos];});

    pop = result;
}

template <typename PHEN_T>
void FilterDominated(emp::vector<PHEN_T> & pop, emp::vector<Axis> & axes, double epsilon = 0) {
    emp::vector<int> to_remove;
    
    for (int i = 0; i < pop.size(); i++) {
        for (int j = 0; j < pop.size(); j++) {
            if (i == j) {
                continue;
            }
            bool losses = false;
            bool wins = false;
            for (Axis & ax : axes) {
                if (pop[i][ax.curr_id] + epsilon < pop[j][ax.curr_id]) {
                    losses = true;
                } else if (pop[i][ax.curr_id] - epsilon > pop[j][ax.curr_id]) {
                    wins = true;
                    break;
                } 
            }
            if (losses && !wins) {
                // std::cout << "Removing " << i << std::endl;
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
    // if (to_remove.size() > 0) {
    // std::cout << "Final pop " << pop << std::endl;
    // }

}


template <typename PHEN_T>
void PruneAxes(emp::vector<Axis> & axes, emp::vector<PHEN_T> & pop, double epsilon = 0) {
    
    emp::vector<int> to_remove;

    for (int i = 0; i < axes.size(); i++) {
        Axis ax = axes[i];
        // std::cout << "Evaluating " << ax << std::endl;
        bool include = false;
        double best = pop[0][ax.curr_id];
        double lowest = pop[0][ax.curr_id];
        for (PHEN_T & org : pop) {
            if (org[ax.curr_id] > best) {
                // std::cout << org[ax] << " > "<< best << std::endl;
                if (org[ax.curr_id] - epsilon > lowest) {
                    // std::cout <<" Including  "<< ax << std::endl;
                    include = true;
                    break;
                }
                best = org[ax.curr_id];
            } else if (org[ax.curr_id] < lowest) {
                // std::cout << org[ax] << " < "<< lowest << std::endl;
                if (org[ax.curr_id] + epsilon < best) {
                    // std::cout <<" Including  "<< ax << std::endl;
                    include = true;
                    break;
                }
                lowest = org[ax.curr_id];
            }
        }

        if (!include) {
            to_remove.push_back(i);
        }
    }

    for (emp::vector<int>::reverse_iterator i = to_remove.rbegin(); i != to_remove.rend(); ++i ) {
        axes.erase(axes.begin()+(*i));
    }
}

template <typename PHEN_T>
void DeDuplicateAxes(emp::vector<Axis> & axes, emp::vector<PHEN_T> & pop) {
    emp::vector<int> to_remove;
    for (int ax = 0; ax < axes.size(); ax++) {
        if (emp::Has(to_remove, ax)) {
            continue;
        }
        for (int ax2 = 0; ax2 < axes.size(); ax2++) {
            if (emp::Has(to_remove, ax2)) {
                continue;
            }
            if (ax == ax2) {
                continue;
            }
            bool same = true;
            for (PHEN_T & org : pop) {
                if (org[axes[ax].curr_id] != org[axes[ax2].curr_id]) {
                    same = false;
                    break;
                }
            }
            if (same) {
                to_remove.push_back(ax2);
                axes[ax].dups += axes[ax2].dups;
                // dups[ax] += dups[ax2];
            }
        }
    }
    std::sort(to_remove.begin(), to_remove.end());
    for (emp::vector<int>::reverse_iterator i = to_remove.rbegin(); i != to_remove.rend(); ++i ) {
        axes.erase(axes.begin()+(*i));
    }

}

template <typename PHEN_T>
void HandleTwoOrgs(std::unordered_map<int, double> & fit_map, emp::vector<PHEN_T> & winners, emp::vector<Axis> & axes, emp::vector<int> & perm_levels, double epsilon = 0, double multiplier=1.0) {
    double wins = 0;
    double losses = 0;
    for (Axis & ax : axes) {
        if (winners[0][ax.curr_id] > winners[1][ax.curr_id] + epsilon) {
            wins += ax.dups;
        } else if (winners[1][ax.curr_id] > winners[0][ax.curr_id] + epsilon) {
            losses += ax.dups;
        }
    }
    if (wins > 0 || losses > 0) {
        fit_map[winners[0].id] += (1.0/winners[0].dups) * (wins/(wins+losses)) * (multiplier/VectorProduct(perm_levels));//(double)emp::Factorial(axes.size() - 1);
        fit_map[winners[1].id] += (1.0/winners[1].dups) * (losses/(wins+losses)) * (multiplier/VectorProduct(perm_levels));//(double)emp::Factorial(axes.size() - 1);
    } else {
        fit_map[winners[0].id] += (1.0/winners[0].dups) * .5 * (multiplier/VectorProduct(perm_levels));//(double)emp::Factorial(axes.size() - 1);
        fit_map[winners[1].id] += (1.0/winners[1].dups) * .5 * (multiplier/VectorProduct(perm_levels));//(double)emp::Factorial(axes.size() - 1);                
    }
}

template <typename PHEN_T>
double HandleTwoOrgsIndividual(emp::vector<PHEN_T> & winners, emp::vector<Axis> & axes, emp::vector<int> & perm_levels, double epsilon = 0) {
    double wins = 0;
    double losses = 0;

    for (Axis & ax : axes) {
        if (winners[0][ax.curr_id] > winners[1][ax.curr_id] + epsilon) {
            wins++;
        } else if (winners[1][ax.curr_id] > winners[0][ax.curr_id] + epsilon) {
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
double TraverseDecisionTreeIndividual(emp::vector<PHEN_T> & pop, emp::vector<Axis> axes, PHEN_T individual, emp::vector<int> perm_levels, double epsilon = 0) {
    // std::cout << "begining round of recursion " << axes.size() << emp::to_string(pop) << emp::to_string(axes) << std::endl;

    // std::cout << emp::to_string(pop) << emp::to_string(axes) << std::endl;
 
    emp_assert(pop.size() > 0, axes, perm_levels);

    // There's only one fitness criterion left, so it wins
    if (axes.size() == 1) {
        emp::vector<PHEN_T> winners = FindHighest(pop, axes[0].curr_id, epsilon);
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
    for (Axis & ax : FindWinningAxes(pop, axes, individual, epsilon)) {
        emp::vector<PHEN_T> winners = FindHighest(pop, ax.curr_id, epsilon);
        emp::vector<Axis> next_axes = axes;
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

template <typename T>
void RankTransform(emp::vector<T> & arr) {
    std::map<T, int> m;
    int count=0;

    for(T i : arr) {
        m[i]++;
    }

    for(auto el : m) {
        m[el.first] = count;
        count++;
    }

    for(int i=0;i<arr.size();i++) {
        arr[i] = m[arr[i]];
    }
}


template <typename PHEN_T>
void Rescale(emp::vector<PHEN_T> & pop, emp::vector<Axis> & axes) {
    emp::vector<typename PHEN_T::value_type> rankOrder(pop.size(), 0);
    emp::vector<typename PHEN_T::value_type> column;

    for (Axis & col : axes) {
        column.clear();
        for (size_t row = 0; row < pop.size(); row++) {
            column.push_back(pop[row][col.curr_id]);
        }

        RankTransform(column);

        for (size_t row = 0; row < pop.size(); row++) {
            pop[row][col.curr_id] = column[row];
        }

    }
}

template <typename PHEN_T>
bool CheckWorstCase(emp::vector<PHEN_T> & pop, emp::vector<Axis> & axes, emp::vector<int> & n_losses) {
    // std::cout << "Starting" << std::endl;
    // bool binary = true;
    for (Axis & ax : axes) {
        int n_change = 0;
        int n_tie = 0;
        int n_lose = 0;
        double best = std::numeric_limits<double>::lowest();
        int index = 0;
        for (PHEN_T & p : pop) {
            // if (p[ax] != 0 && p[ax] != 1) {
            //     binary = false;
            // }
            if (p[ax.curr_id] > best) {
                if (n_tie > 1) {
                    return false;
                }
                if (n_tie == 1) {
                    // if (n_losses[0] != 0) {
                    //     return false;
                    // }
                    n_losses[0] += ax.dups;
                }
                best = p[ax.curr_id];
                n_change++;
                n_lose = n_tie;
                n_tie = 1;
            } else if (almost_equal((double) p[ax.curr_id], best, 5)) {
                n_tie++;
            } else {
                // if (n_losses[index] != 0) {
                //     return false;
                // }
                n_losses[index] += ax.dups;
                n_lose++;
            }
            if (n_change > 2) {
                return false;
            }
            if (n_lose > 1) {
                return false;
            }
            index++;
        }
    }
    std::cout << "TRUE" << std::endl;
    return true;
}

long long unsigned int SwitchFactorial(int i) {
    switch (i)
    {
    case 0:
        return 1;
    case 1:
        return 1;

    case 2: 
        return 2;

    case 3:
        return 6;

    case 4:
        return 24;

    case 5: 
        return 120;
    
    case 6:
        return 720;

    case 7:
        return 5040;

    case 8:
        return 40320;

    case 9:
        return 362880;

    case 10:
        return 3628800;

    case 11:
        return 39916800;

    case 12:
        return 479001600;

    case 13:
        return 6227020800;

    case 14:
        return 87178291200;

    case 15:
        return 1307674368000;

    case 16:
        return 20922789888000;

    case 17:
        return 355687428096000;

    case 18:
        return 6402373705728000;

    case 19:
        return 121645100408832000;

    case 20:
        return 2432902008176640000;
  
    default:
        emp_assert(false, i);
        return -1;
    }
}

constexpr long long unsigned int Factorial(int i) {
    emp_assert(i < 21, i);
    long long unsigned int result = 1;
    while (i > 0) {
      result *= i;
      i--;
    }
    return result;
}

constexpr long long unsigned int FactorialDiff(int i, int j) {
    if (i==0){
        return 1;
    }
    long long unsigned int result = 1;
    for (int curr = i; curr > j; curr--) {
        result *= curr;
    }
    return result;
  }


long long unsigned int DoCombinatorics(int A, int B, int C, int E, int F, int G) {
    long long unsigned int result = 0;
    int total = A + B + C + E + F + G;
    for (int i = 0; i < A+E; i++) {
        result += A * FactorialDiff(A+E-1, A+E-1-i) * (B+C) * Factorial(total-2-i);
    }
    for (int i = 0; i < C+F; i++) {
        result += C * FactorialDiff(C+F-1, C+F-1-i) * (B+A) * Factorial(total-2-i);
    }

    result += B * Factorial(total-1);
    return result;
}


template <typename PHEN_T>
void HandleThreeOrgs(std::unordered_map<int, double> & fit_map, emp::vector<PHEN_T> & winners, emp::vector<Axis> & axes, emp::vector<int> & perm_levels, double multiplier=1.0) {
    using inner_type = typename PHEN_T::value_type;
    using vals_type = emp::vector<inner_type>;
    // std::cout << winners << std::endl;
    std::map<vals_type, int> counts;
    
    // counts[{0,0,0}] = 0; // Can't actually happen
    // counts[{0,0,1}] = 0;
    // counts[{0,1,0}] = 0;
    // counts[{1,0,0}] = 0;
    // counts[{0,1,1}] = 0;
    // counts[{1,0,1}] = 0;
    // counts[{1,1,0}] = 0;
    // counts[{1,1,1}] = 0; // Can't actually happen

    // Can ignore everything that is just 2s and
    // 1s or 2s and 0s because those can be mapped
    // to combos of 0s and 1s

    // counts[{0,1,2}] = 0;
    // counts[{0,2,1}] = 0;
    // counts[{1,0,2}] = 0;
    // counts[{1,2,0}] = 0;
    // counts[{2,1,0}] = 0;
    // counts[{2,0,1}] = 0;

    int total = 0;

    for (Axis & ax : axes) {
        vals_type vals;
        for (PHEN_T & p : winners) {
            vals.push_back(p[ax.curr_id]);
        }

        // if (emp::Has<inner_type>(vals, 0) && emp::Has<inner_type>(vals, 2) && !emp::Has<inner_type>(vals, 1)) {
        //     std::replace_if(vals.begin(), vals.end(), [](inner_type & i){return i==2;}, 1);
        // }
        // if (emp::Has<inner_type>(vals, 1) && emp::Has<inner_type>(vals, 2) && !emp::Has<inner_type>(vals, 0)) {
        //     std::replace_if(vals.begin(), vals.end(), [](inner_type & i){return i==1;}, 0);
        //     std::replace_if(vals.begin(), vals.end(), [](inner_type & i){return i==2;}, 1);
        // }

        counts[vals] += ax.dups;
        // std::cout << vals << " " << counts[vals] << std::endl;
        total += ax.dups;
    }

    // for (auto el : counts) {
    //     std::cout << el.first << " " << el.second << std::endl;
    // }
    
    long long unsigned int org0 = DoCombinatorics(counts[{1, 0, 1}],
                                counts[{1, 0, 0}],
                                counts[{1, 1, 0}],
                                counts[{0, 1, 0}],
                                counts[{0, 0, 1}],
                                counts[{0, 1, 1}] 
    
    );
    long long unsigned int org1 = DoCombinatorics(counts[{0, 1, 1}],
                               counts[{0, 1, 0}],
                               counts[{1, 1, 0}],
                               counts[{1, 0, 0}],
                               counts[{0, 0, 1}], 
                               counts[{1, 0, 1}]
    );
    long long unsigned int org2 = DoCombinatorics(counts[{0, 1, 1}],
                               counts[{0, 0, 1}], 
                               counts[{1, 0, 1}],
                               counts[{1, 0, 0}],
                               counts[{0, 1, 0}],
                               counts[{1, 1, 0}]
    );

    emp_assert(counts.size() <= 6, counts.size());
    int empty_count = counts[{0,0,0}];
    emp_assert(empty_count == 0, empty_count);
    empty_count = counts[{1,1,1}];
    emp_assert(empty_count == 0, empty_count);
    emp_assert((org0 + org1 + org2) == Factorial(total), org0, org1, org2, org0 + org1 + org2, Factorial(total));
    // fit_map[winners[0].id]+= (multiplier/winners[0].dups) * 1.0/emp::from_string<int>((big_fact(total)/org0).GetString()) * 1.0/(1 * VectorProduct(perm_levels));
    // fit_map[winners[1].id]+= (multiplier/winners[1].dups) * 1.0/emp::from_string<int>((big_fact(total)/org1).GetString()) * 1.0/(1 * VectorProduct(perm_levels));        
    // fit_map[winners[2].id]+= (multiplier/winners[2].dups) * 1.0/emp::from_string<int>((big_fact(total)/org2).GetString()) * 1.0/(1 * VectorProduct(perm_levels));
    fit_map[winners[0].id]+= (multiplier/winners[0].dups) * ((double)org0)/Factorial(total) * 1.0/(1 * VectorProduct(perm_levels));
    fit_map[winners[1].id]+= (multiplier/winners[1].dups) * ((double)org1)/Factorial(total) * 1.0/(1 * VectorProduct(perm_levels));        
    fit_map[winners[2].id]+= (multiplier/winners[2].dups) * ((double)org2)/Factorial(total) * 1.0/(1 * VectorProduct(perm_levels));

    // std::cout << org0 << " " << org1 <<  " " << org2 << " " << big_fact(total) << std::endl;
}

template <typename PHEN_T>
void TraverseDecisionTree(std::unordered_map<int, double> & fit_map, emp::vector<PHEN_T> & pop, emp::vector<Axis> & axes, emp::vector<int> perm_levels, double epsilon = 0, double multiplier = 1.0, bool binary = false) {
    // std::cout << "begining round of recursion " << axes.size() << emp::to_string(pop) << emp::to_string(axes) << std::endl;

    // std::cout << emp::to_string(pop) << emp::to_string(axes) << std::endl;
 
    emp_assert(pop.size() > 0, pop,  axes, perm_levels);
    emp_assert(axes.size() > 0, axes, perm_levels);

    // There's only one fitness criterion left, so it wins
    if (axes.size() == 1) {
        emp::vector<PHEN_T> winners = FindHighest(pop, axes[0].curr_id, epsilon);
        // std::cout << "Winners: " << emp::to_string(winners) << std::endl;
        for (PHEN_T & org : winners) {
            fit_map[org.id]+= (1.0/org.dups) * multiplier*axes[0].dups/((double)winners.size()*VectorProduct(perm_levels));
        }
        return;
    }

    // There's only one thing in the population, so it wins
    if (pop.size() == 1) {
        // TODO: Figure out if we actually need this
        // emp::vector<PHEN_T> winners = FindHighest(pop, axes[0], epsilon);
        // std::cout << "Winners: " << emp::to_string(winners) << std::endl;
        // for (PHEN_T & org : pop) {
            // std::cout << org.curr_phen << " " << org.dups << " " << (1.0/org.dups) * multiplier/VectorProduct(perm_levels) << std::endl;
        fit_map[pop[0].id]+= (1.0/pop[0].dups) * multiplier/VectorProduct(perm_levels);
        // }
        return;
    }

    PruneAxes(axes, pop, epsilon);
    if (axes.size() > pop.size()) {
        Rescale(pop, axes);
    }    
    DeDuplicateAxes(axes, pop);

    int real_axes = 0;
    for (Axis & ax : axes){
        real_axes += ax.dups;
    }

    if (axes.size() == 0) {
        for (PHEN_T & org : pop) {
            fit_map[org.id]+= (1.0/org.dups) * multiplier/((double)pop.size()*VectorProduct(perm_levels));
        }        
        return;
    }

    perm_levels.push_back(real_axes);

    // Check for only one axis again now that we've pruned the axes
    if (axes.size() == 1) {
        emp::vector<PHEN_T> winners = FindHighest(pop, axes[0].curr_id, epsilon);
        // std::cout << "Winners: " << emp::to_string(winners) << std::endl;
        for (PHEN_T & org : winners) {
            fit_map[org.id]+=(1.0/org.dups) * multiplier*real_axes*axes[0].dups/((double)winners.size()*VectorProduct(perm_levels));
        }
        return;
    }

    // std::cout << "post processing: " << axes.size() << emp::to_string(pop) << emp::to_string(axes) << std::endl;
    // emp::vector<int> n_losses(pop.size(), 0);
    // if (CheckWorstCase(pop, axes, n_losses)) {
    //     // std::cout << n_losses << std::endl;
    //     // emp_assert(emp::Sum(n_losses) == axes.size(), n_losses, axes, pop);
    //     for (int i = 0; i < pop.size(); i++) {
    //         // int win_axes = real_axes - n_losses[i];
    //         // int won_dups = real_axes - axes.size() - n_losses[i];
    //         // int wins = emp::Factorial(win_axes);
    //         // for (int j = 0; j < won_dups; j++) {
    //         //     win_axes
    //         // }
    //         fit_map[pop[i].id] += (1.0/pop[i].dups) * multiplier/((double)pop.size()*VectorProduct(perm_levels));
    //     }
    //     return;
    // }

    for (int i = 0; i < axes.size(); i++) {
        Axis ax = axes[i];
        // if (perm_levels.size() == 1) {
        //     std::cout << "Axis: " << ax << " out of " << emp::to_string(axes) << std::endl;
        // }
        emp::vector<PHEN_T> winners = FindHighest(pop, ax.curr_id, epsilon);
        emp::vector<Axis> next_axes = axes;
        next_axes.erase(next_axes.begin()+i);

        if (epsilon) {
            FilterImpossible(winners, next_axes, epsilon);
        } else {
            FilterDominated(winners, next_axes, epsilon);
            if (winners.size() > 10) {
                Rescale(winners, axes);
                // winners = RemoveDuplicatesIndexes(winners, next_axes);
            }
        }
        
        double new_multiplier = multiplier * ax.dups;
        int num_axes = 0;
        for (Axis & ax_i : next_axes) {
            num_axes += ax_i.dups;
        }

        if (winners.size() == 1) { // Not a tie
            // std::cout << "1 winner: " << emp::to_string(winners[0]) << " Controls " << (double)emp::Factorial(axes.size() - 1)<< std::endl;
            fit_map[winners[0].id] += (1.0/winners[0].dups) *new_multiplier/VectorProduct(perm_levels);//(double)emp::Factorial(axes.size() - 1);
        } else if (winners.size() == 2) { // optimization
            HandleTwoOrgs(fit_map, winners, next_axes, perm_levels, epsilon, new_multiplier);
        // } else if (winners.size() < 100 && next_axes.size() > 2) { // optimization
        //     HandleThreeOrgs(fit_map, winners, next_axes, perm_levels, epsilon);
        } else if (binary && winners.size() == 3 && epsilon == 0 && num_axes < 18) {
            PruneAxes(next_axes, winners, epsilon);
            Rescale(winners, next_axes);
            // DeDuplicateAxes(next_axes, winners);
            // winners = RemoveDuplicatesIndexes(winners, next_axes);

            // std::cout << "Axes " << next_axes << std::endl;
            HandleThreeOrgs(fit_map, winners, next_axes, perm_levels, new_multiplier);
        } else { // tie
            TraverseDecisionTree(fit_map, winners, next_axes, perm_levels, epsilon, new_multiplier, binary);
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
emp::vector<Individual<PHEN_T>> SetupPop(emp::vector<PHEN_T> & pop, std::map<PHEN_T, int> & seen) {
    emp::vector<Individual<PHEN_T>> new_pop;
    for (PHEN_T p : pop) {
        if (emp::Has(seen, p)) {
            new_pop[seen[p]].dups++;
        } else {
            seen[p] = new_pop.size();
            new_pop.emplace_back(p, new_pop.size());
        }
    }
    return new_pop;
}


template <typename PHEN_T>
emp::vector<double> LexicaseFitness(emp::vector<PHEN_T> & pop, double epsilon = 0, bool binary=false) {

    if (pop.size() == 0) {
        return emp::vector<double>();
    }

    std::map<PHEN_T, int> phen_to_id;
    emp::vector<Individual<PHEN_T>> de_dup_pop = SetupPop(pop, phen_to_id);

    std::unordered_map<int, double> fit_map;
    size_t n_funs = pop[0].size();

    for (Individual<PHEN_T> & org : de_dup_pop) {
        fit_map[org.id] = 0.0;
    }

    emp::vector<Axis> axes;
    for (int i = 0; i < n_funs; i++) {
        axes.emplace_back(i);
    }
    // emp::vector<PHEN_T> de_dup_pop = emp::RemoveDuplicates(pop);
    if (epsilon) {
        FilterImpossible(de_dup_pop, axes, epsilon);
    } else {
        FilterDominated(de_dup_pop, axes, epsilon);
    }    
    if (binary && de_dup_pop.size() == 3 && epsilon == 0) {
        PruneAxes(axes, de_dup_pop, epsilon);
        Rescale(de_dup_pop, axes);
        DeDuplicateAxes(axes, de_dup_pop);
        // winners = RemoveDuplicatesIndexes(winners, next_axes);

        // std::cout << "Axes " << next_axes << std::endl;
        emp::vector<int> perm_levels(0);
        HandleThreeOrgs(fit_map, de_dup_pop, axes, perm_levels, 1.0);
    } else {
        TraverseDecisionTree(fit_map, de_dup_pop, axes, {}, epsilon, 1.0, binary);
    }


    // for (PHEN_T & org : de_dup_pop) {
    //     fit_map[org] /= emp::Count(pop, org);
    // }

    emp::vector<double> result;
    for (PHEN_T & org : pop) {
        // std::cout << org << " " << phen_to_id[org] << " " << fit_map[phen_to_id[org]] << std::endl;
        result.push_back(fit_map[phen_to_id[org]]);
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
    emp::vector<Axis> axes;
    for (int i = 0; i < n_funs; i++) {
        axes.emplace_back(i);
    }
    return TraverseDecisionTreeIndividual(de_dup_pop, axes, phen, {}, epsilon)/n_dups;
}


void TournamentHelper(emp::vector<double> & fit_map, int t_size = 2){

    emp::vector<double> base_fit_map = fit_map;
    double pop_size = base_fit_map.size();

    for (size_t i = 0; i < fit_map.size(); i++) {
        double less = 0.0;
        double equal = 0.0; 

        for (auto & org2 : base_fit_map) {
            if (almost_equal(org2, base_fit_map[i], 10)) {
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


template <typename PHEN_T>
emp::WeightedGraph CalcCompetition(emp::vector<PHEN_T> pop, 
                        std::function<emp::vector<double>(emp::vector<PHEN_T>&)> fit_fun) {

    emp::WeightedGraph effects(pop.size());


    emp::vector<double> fitnesses = fit_fun(pop);

    for (size_t i = 0; i < pop.size(); i++) {
        effects.SetLabel(i, emp::to_string(pop[i]));
        // std::cout << effects.GetLabel(i) << std::endl;

        emp::vector<PHEN_T> curr = pop;
        for (size_t ax = 0; ax < curr[i].size(); ax++) {
            curr[i][ax] = 0; // Replace org with null org so pop size stays same
        }

        emp::vector<double> new_fits = fit_fun(curr);
        for (size_t j = 0; j < pop.size(); j++ ) {
            if (i == j) {continue;}

            // In terms of floating-point imprecision issues, it's much better
            // to check for equality than doing the subtraction and checking
            // for the result to be 0
            if (!almost_equal(fitnesses[j], new_fits[j], 10)) {
                effects.AddEdge(i, j, fitnesses[j] - new_fits[j]);
            }
            // std::cout << effect << std::endl;
        }
    }

    return effects;
}



struct Node {
    std::set<Axis> axes;
    std::map<int, long long unsigned int> seq_counts;
    std::set<emp::Ptr<Node>> children;
    int free_axes = 0;
    int real_axes = 0;
    Node(std::set<Axis> ax) {
        axes = ax;
    }
    Node() {}
};

template <typename PHEN_T>
int SetupBinaryPop(emp::vector<PHEN_T> & pop, std::map<PHEN_T, int> & counts, int ind) {
    emp::vector<PHEN_T> new_pop;
    int new_ind = -1;
    for (int i = 0; i < pop.size(); i++) {
        if (emp::Has(counts, pop[i])) {
            counts[pop[i]]++;
        } else {
            if (pop[i] == pop[ind]) {
                new_ind = new_pop.size();
            }
            counts[pop[i]] = 1;
            new_pop.push_back(pop[i]);
        }
    }
    pop = new_pop;
    return new_ind;
}

template <typename PHEN_T>
double SolveBinary(emp::vector<PHEN_T> & pop, int ind) {
    using inner_type = typename PHEN_T::value_type;
    using vals_type = emp::vector<inner_type>;

    std::map<PHEN_T, int> counts;
    ind = SetupBinaryPop(pop, counts, ind);
    int dup_solutions = counts[pop[ind]];

    std::map<std::set<Axis>, Node> nodes;
    emp::vector<Axis> all_axes;
    for (int i = 0; i < pop[ind].size(); i++) {
        all_axes.emplace_back(i);
    }
    PruneAxes(all_axes, pop);
    int n_real_axes = all_axes.size();
    DeDuplicateAxes(all_axes, pop);

    for (Axis & ax1 : all_axes) {
        if (!pop[ind][ax1.curr_id]) {
            continue;
        }
        for (Axis & ax2 : all_axes) {
            if (ax1 == ax2) {
                continue;
            }
            bool subset = true;
            for (int i = 0; i < pop.size(); i++) {
                if (i == ind) {
                    continue;
                }
                if (pop[i][ax1.curr_id] && !pop[i][ax2.curr_id]) {
                    subset = false;
                }
            }
            if (subset) {
                ax1.sub_axes.push_back(ax2.curr_id);
                ax1.bonus_dups += ax2.dups;
            }

        }
    }

    auto cmp = [](emp::Ptr<Node> left, emp::Ptr<Node> right) { return left->axes.size() > right->axes.size(); };
    std::priority_queue<emp::Ptr<Node>, std::vector<emp::Ptr<Node>>, decltype(cmp)> pq(cmp);

    // emp::vector<emp::Ptr<Node>> start_nodes(0);
    for (Axis & ax : all_axes) {
        if (pop[ind][ax.curr_id]) {
            nodes[{ax}] = Node({ax});
            emp::Ptr<Node> curr_node = &nodes[{ax}];
            curr_node->real_axes = 1;
            pq.push(curr_node);
            curr_node->axes.insert(ax.sub_axes.begin(), ax.sub_axes.end());
            for (int i = 1; i <= ax.dups + ax.bonus_dups; i++) {
                int avail = ax.dups + ax.bonus_dups - 1;
                curr_node->seq_counts[i] = ax.dups * FactorialDiff(avail, avail - i + 1);
            }
            curr_node->free_axes = ax.dups + ax.bonus_dups - 1;
        }
    }


    // emp::vector<emp::Ptr<Node>> curr_nodes = start_nodes;
    // emp::vector<emp::Ptr<Node>> next_nodes;
    // for (int depth = 2; depth <= all_axes.size(); depth++) {
    //     for (emp::Ptr<Node> n : curr_nodes ) {
    while(pq.size()) {
        std::cout << pq.size() << std::endl;
        emp::Ptr<Node> n = pq.top();
        pq.pop();
        for (Axis & ax : all_axes) {
            if (emp::Has(n->axes, ax)) {
                continue;
            }
            // bool can_use = true;
            if (!pop[ind][ax.curr_id]) {
                // for (int j = 0; j < pop.size(); j++) {
                //     if (pop[j][ax.curr_id] && std::all_of(n->axes.begin(), n->axes.end(), [j, &pop](const Axis & x){return pop[j][x.curr_id];})) {
                //         can_use = false;
                //         break;
                //     }
                // }
                // Identify members of current set that have
                // 1s for individual
                // Try choosing all of those, then all of the others,
                // then new axis. See if that actually works.
                // if (can_use) {
                emp::vector<Axis> ones;
                emp::vector<Axis> zeros;

                for (Axis ax2 : n->axes) {
                    if (pop[ind][ax2.curr_id]) {
                        ones.push_back(ax2);
                    } else {
                        zeros.push_back(ax2);
                    }
                }
                emp::vector<PHEN_T> test_pop = pop;
                for (Axis & ax2 : ones) {
                    test_pop = FindHighest(test_pop, ax2.curr_id);
                }
                for (Axis & ax2 : zeros) {
                    test_pop = FindHighest(test_pop, ax2.curr_id);
                }
                test_pop = FindHighest(test_pop, ax.curr_id);
                if (!emp::Has(test_pop, pop[ind])) {
                    // can_use = false;
                    continue;
                }

            }


            std::set<Axis> temp_set = n->axes;
            temp_set.insert(ax);
            temp_set.insert(ax.sub_axes.begin(), ax.sub_axes.end());

            if (!emp::Has(nodes, temp_set)) {
                nodes[temp_set] = Node(temp_set);
                pq.push(&nodes[temp_set]);
                nodes[temp_set].free_axes = n->free_axes + ax.dups + ax.bonus_dups - 1;
            }
            emp::Ptr<Node> curr_node = &nodes[temp_set];
            n->children.insert(curr_node);
            // int max_val = std::max_element(n->seq_counts.begin(), n->seq_counts.end(), [](auto & p1, auto & p2){return p1.first < p2.first;})->first;
            int avail = ax.dups - 1 + ax.bonus_dups + n->free_axes;
            for (int i = 1; i < n->axes.size() + avail; i++) {
                // depth + 1 (dups options)
                // depth + 2 (dups * free axes) - free_axes is max_val - depth
                // curr_node->seq_counts[i] += n->seq_counts[i-1] * ax.dups;
                
                for (int j = 1; j < i; j++){
                    curr_node->seq_counts[i] += n->seq_counts[i-j] * ax.dups * FactorialDiff(avail, avail - j + 1); 
                }
            }
            // curr_node->free_axes += ax.dups - 1;
        }
        // }
        // for (emp::Ptr<Node> n : curr_nodes) {
        //     nodes.erase(n->axes);
        // }
        // curr_nodes = next_nodes;
        // next_nodes.clear();
    }
    Node & answer = nodes[std::set(all_axes.begin(), all_axes.end())];
    // emp_assert(curr_nodes.size() == 1, curr_nodes.size());
    return ((double)answer.seq_counts[n_real_axes]/dup_solutions)/(double)SwitchFactorial(n_real_axes);
}



#endif
