import ec_ecology_toolbox.selection_probabilities as eco
import random
import pytest


def test_lex_prob():
    assert eco.LexicaseFitness([[1, 2, 3]], 1) == [1]
    assert eco.LexicaseFitness([[1, 2, 3], [2, 1, 4]], 1) == [.5, .5]
    result = eco.LexicaseFitness([[1, 2, 3], [2, 1, 4]])
    assert result[0] == pytest.approx(.3333333333)
    assert result[1] == pytest.approx(.6666666667)
    assert eco.LexicaseFitness([]) == []
    assert eco.LexicaseFitness([[]]) == [1]


@pytest.mark.parametrize("seed, n, m",
                         [(1, 10, 10),
                          (2, 10, 20),
                          (3, 10, 30),
                          (4, 20, 10),
                          (5, 20, 20),
                          (6, 3, 20),                          
                          (7, 20, 3),                          
                          (8, 8, 20),                          
                          (9, 20, 8),                            
                        #   (6, 20, 30),
                        #   (7, 30, 30),
                        #   (8, 40, 40)
                          ])
def test_benchmark_lex_prob(benchmark, seed, n, m):
    random.seed(seed)
    benchmark(eco.LexicaseFitness,
              [[random.randint(0, 1) for i in range(n)] for j in range(m)], 0)

# @pytest.mark.parametrize("seed, n, m",
#                          [(1, 10, 10),
#                           (2, 10, 20),
#                           (3, 10, 30),
#                           (4, 20, 10),
#                           (5, 20, 20),
#                           (6, 3, 20),                          
#                           (7, 20, 3),                          
#                           (8, 8, 20),                          
#                           (9, 20, 8),                          
#                         #   (6, 20, 30),
#                         #   (7, 30, 30),
#                         #   (8, 40, 40)
#                           ])
# def test_benchmark_binary_lex_prob(benchmark, seed, n, m):
#     random.seed(seed)
#     benchmark(eco.LexicaseFitnessIndividualBinary,
#               [[random.randint(0, 1) for i in range(n)] for j in range(m)], 0)


def test_lex_prob_individual():
    assert eco.LexicaseFitnessIndividual([[1, 2, 3]], 0, 1) == 1
    assert eco.LexicaseFitnessIndividual([[1, 2, 3], [2, 1, 4]], 1, 1.1) == .5
    assert eco.LexicaseFitnessIndividual([[1, 2, 3], [2, 1, 4]], 0) == pytest.approx(.3333333333)
    assert eco.LexicaseFitnessIndividual([[1, 2, 3], [2, 1, 4]], 1) == pytest.approx(.6666666667)


def test_sharing_prob():
    result = eco.SharingFitness([[1, 2, 3], [1, 2, 3], [3, 2, 1]])
    assert result[2] == pytest.approx(5/9)
    assert result[0] == pytest.approx(2/9)
    assert result[1] == pytest.approx(2/9)

    result = eco.SharingFitness([[1, 2, 3], [3, 2, 1]])
    assert result[0] == pytest.approx(.5)
    assert result[1] == pytest.approx(.5)

    result = eco.SharingFitness([[1, 2, 3], [1, 2, 3]])
    assert result[0] == pytest.approx(.5)
    assert result[1] == pytest.approx(.5)

    result = eco.SharingFitness([[1, 2, 3], [2, 1, 3], [3, 2, 1]], sigma_share=1)
    assert result[2] == pytest.approx(3/9)
    assert result[0] == pytest.approx(3/9)
    assert result[1] == pytest.approx(3/9)

    result = eco.SharingFitness([[1, 2, 3], [2, 1, 3], [3, 2, 1]], sigma_share=2)
    assert result[2] == pytest.approx(5/9)
    assert result[0] == pytest.approx(2/9)
    assert result[1] == pytest.approx(2/9)

    result = eco.SharingFitness([[1, 2, 3], [2, 1, 3], [3, 2, 1]], t_size=1)
    assert result[2] == pytest.approx(3/9)
    assert result[0] == pytest.approx(3/9)
    assert result[1] == pytest.approx(3/9)

    result = eco.SharingFitness([[1, 2, 3], [2, 1, 3], [3, 2, 1]], sigma_share=2, t_size = 3)
    assert result[2] == pytest.approx(1 - (2/3)**3)
    assert result[0] == pytest.approx(((2/3)**3)/2)
    assert result[1] == pytest.approx(((2/3)**3)/2)


def test_tournament_prob():
    result = eco.TournamentFitness([[1, 2, 3], [2, 1, 3], [3, 2, 1]], t_size=1)
    assert result[2] == pytest.approx(3/9)
    assert result[0] == pytest.approx(3/9)
    assert result[1] == pytest.approx(3/9)

    result = eco.TournamentFitness([[1, 2, 3], [2, 1, 3], [3, 2, 1]], t_size=2)
    assert result[2] == pytest.approx(3/9)
    assert result[0] == pytest.approx(3/9)
    assert result[1] == pytest.approx(3/9)

    result = eco.TournamentFitness([[1, 2, 3], [2, 1, 3], [3, 2, 1]], t_size=3)
    assert result[2] == pytest.approx(3/9)
    assert result[0] == pytest.approx(3/9)
    assert result[1] == pytest.approx(3/9)

    result = eco.TournamentFitness([[1, 2, 3], [2, 1, 3], [3, 3, 1]], t_size=3)
    assert result[2] == pytest.approx(1 - (2/3)**3)
    assert result[0] == pytest.approx(((2/3)**3)/2)
    assert result[1] == pytest.approx(((2/3)**3)/2)


def test_nk():
    r = eco.Random(-1)
    nk = eco.NKLandscape(8, 2, r)
    print(nk.GetFitness(0, 1))
