import ec_ecology_toolbox as eco
import pytest


def test_lex_prob():
    assert eco.LexicaseFitness([[1, 2, 3]], 1) == [1]
    assert eco.LexicaseFitness([[1, 2, 3], [2, 1, 4]], 1) == [.5, .5]
    result = eco.LexicaseFitness([[1, 2, 3], [2, 1, 4]])
    assert result[0] == pytest.approx(.3333333333)
    assert result[1] == pytest.approx(.6666666667)
    assert eco.LexicaseFitness([]) == []
    assert eco.LexicaseFitness([[]]) == [1]