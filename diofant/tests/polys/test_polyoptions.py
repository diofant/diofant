"""Tests for options manager for :class:`Poly` and public API functions."""

import pytest

from diofant import (CC, EX, FF, GF, QQ, RR, ZZ, ComplexField, GeneratorsError,
                     I, Integer, OptionError, Options, RealField, Symbol, lex,
                     sqrt)
from diofant.abc import x, y, z
from diofant.polys.polyoptions import (All, Auto, BooleanOption, Domain,
                                       Expand, Extension, Field, Formal, Frac,
                                       Gaussian, Gen, Gens, Greedy, Include,
                                       Method, Modulus, OptionType, Order,
                                       Polys, Sort, Split, Strict, Symbols,
                                       Wrt, allowed_flags, set_defaults)


__all__ = ()


def test_Options_clone():
    opt = Options((x, y, z), {'domain': 'ZZ'})

    assert opt.gens == (x, y, z)
    assert opt.domain == ZZ
    assert ('order' in opt) is False
    assert opt.args == {'domain': ZZ}

    # defaults:
    assert opt.flags['all'] is False
    assert opt.flags['include'] is False
    assert opt.options['strict'] is True

    new_opt = opt.clone({'gens': (x, y), 'order': 'lex'})

    assert opt.gens == (x, y, z)
    assert opt.domain == ZZ
    assert ('order' in opt) is False

    assert new_opt.gens == (x, y)
    assert new_opt.domain == ZZ
    assert ('order' in new_opt) is True

    opt.spam = 'eggs'
    assert opt.spam == 'eggs'

    class SpamOpt(BooleanOption, metaclass=OptionType):
        option = 'spam'
        before = ['gens']
        after = ['domain']

    Options.__order__ = None
    pytest.raises(RuntimeError, lambda: Options._init_dependencies_order())
    delattr(Options, 'spam')
    del Options.__options__['spam']
    Options.__order__ = None
    Options._init_dependencies_order()
    Options._init_dependencies_order()  # noop

    pytest.raises(OptionError, lambda: Options((x,), {'gens': (x, y)}))
    pytest.raises(OptionError, lambda: Options((x,), {'spam': 1}))
    pytest.raises(OptionError, lambda: Options((x,), {'field': True,
                                                      'gaussian': True}))
    pytest.raises(OptionError, lambda: Options((x,), {'gen': x}, strict=True))


def test_Expand_preprocess():
    assert Expand.preprocess(False) is False
    assert Expand.preprocess(True) is True

    assert Expand.preprocess(0) is False
    assert Expand.preprocess(1) is True

    pytest.raises(OptionError, lambda: Expand.preprocess(x))


def test_Expand_postprocess():
    opt = {'expand': True}
    Expand.postprocess(opt)

    assert opt == {'expand': True}


def test_Gens_preprocess():
    assert Gens.preprocess((None,)) == ()
    assert Gens.preprocess((x, y, z)) == (x, y, z)

    a = Symbol('a', commutative=False)

    pytest.raises(GeneratorsError, lambda: Gens.preprocess((x, x, y)))
    pytest.raises(GeneratorsError, lambda: Gens.preprocess((x, y, a)))


def test_Gens_postprocess():
    opt = {'gens': (x, y)}
    Gens.postprocess(opt)

    assert opt == {'gens': (x, y)}


def test_Wrt_preprocess():
    assert Wrt.preprocess(x) == ['x']
    assert Wrt.preprocess('') == []
    assert Wrt.preprocess(' ') == []
    assert Wrt.preprocess('x,y') == ['x', 'y']
    assert Wrt.preprocess('x y') == ['x', 'y']
    assert Wrt.preprocess('x, y') == ['x', 'y']
    assert Wrt.preprocess('x , y') == ['x', 'y']
    assert Wrt.preprocess(' x, y') == ['x', 'y']
    assert Wrt.preprocess(' x,  y') == ['x', 'y']
    assert Wrt.preprocess([x, y]) == ['x', 'y']

    pytest.raises(OptionError, lambda: Wrt.preprocess(','))
    pytest.raises(OptionError, lambda: Wrt.preprocess(0))


def test_Wrt_postprocess():
    opt = {'wrt': ['x']}
    Wrt.postprocess(opt)

    assert opt == {'wrt': ['x']}


def test_Sort_preprocess():
    assert Sort.preprocess([x, y, z]) == ['x', 'y', 'z']
    assert Sort.preprocess((x, y, z)) == ['x', 'y', 'z']

    assert Sort.preprocess('x > y > z') == ['x', 'y', 'z']
    assert Sort.preprocess('x>y>z') == ['x', 'y', 'z']

    pytest.raises(OptionError, lambda: Sort.preprocess(0))
    pytest.raises(OptionError, lambda: Sort.preprocess({x, y, z}))


def test_Sort_postprocess():
    opt = {'sort': 'x > y'}
    Sort.postprocess(opt)

    assert opt == {'sort': 'x > y'}


def test_Order_preprocess():
    assert Order.preprocess('lex') == lex


def test_Order_postprocess():
    opt = {'order': True}
    Order.postprocess(opt)

    assert opt == {'order': True}


def test_Field_preprocess():
    assert Field.preprocess(False) is False
    assert Field.preprocess(True) is True

    assert Field.preprocess(0) is False
    assert Field.preprocess(1) is True

    pytest.raises(OptionError, lambda: Field.preprocess(x))


def test_Field_postprocess():
    opt = {'field': True}
    Field.postprocess(opt)

    assert opt == {'field': True}


def test_Greedy_preprocess():
    assert Greedy.preprocess(False) is False
    assert Greedy.preprocess(True) is True

    assert Greedy.preprocess(0) is False
    assert Greedy.preprocess(1) is True

    pytest.raises(OptionError, lambda: Greedy.preprocess(x))


def test_Greedy_postprocess():
    opt = {'greedy': True}
    Greedy.postprocess(opt)

    assert opt == {'greedy': True}


def test_Domain_preprocess():
    assert Domain.preprocess(ZZ) == ZZ
    assert Domain.preprocess(QQ) == QQ
    assert Domain.preprocess(EX) == EX
    assert Domain.preprocess(FF(2)) == FF(2)
    assert Domain.preprocess(ZZ.inject(x, y)) == ZZ.inject(x, y)

    assert Domain.preprocess('Z') == ZZ
    assert Domain.preprocess('Q') == QQ

    assert Domain.preprocess('ZZ') == ZZ
    assert Domain.preprocess('QQ') == QQ

    assert Domain.preprocess('EX') == EX

    assert Domain.preprocess('FF(23)') == FF(23)
    assert Domain.preprocess('GF(23)') == GF(23)

    pytest.raises(OptionError, lambda: Domain.preprocess('Z[]'))

    assert Domain.preprocess('Z[x]') == ZZ.inject(x)
    assert Domain.preprocess('Q[x]') == QQ.inject(x)

    assert Domain.preprocess('ZZ[x]') == ZZ.inject(x)
    assert Domain.preprocess('QQ[x]') == QQ.inject(x)

    assert Domain.preprocess('Z[x,y]') == ZZ.inject(x, y)
    assert Domain.preprocess('Q[x,y]') == QQ.inject(x, y)

    assert Domain.preprocess('ZZ[x,y]') == ZZ.inject(x, y)
    assert Domain.preprocess('QQ[x,y]') == QQ.inject(x, y)

    pytest.raises(OptionError, lambda: Domain.preprocess('Z()'))

    assert Domain.preprocess('Z(x)') == ZZ.inject(x).field
    assert Domain.preprocess('Q(x)') == QQ.inject(x).field

    assert Domain.preprocess('ZZ(x)') == ZZ.inject(x).field
    assert Domain.preprocess('QQ(x)') == QQ.inject(x).field

    assert Domain.preprocess('Z(x,y)') == ZZ.inject(x, y).field
    assert Domain.preprocess('Q(x,y)') == QQ.inject(x, y).field

    assert Domain.preprocess('ZZ(x,y)') == ZZ.inject(x, y).field
    assert Domain.preprocess('QQ(x,y)') == QQ.inject(x, y).field

    assert Domain.preprocess('Q<I>') == QQ.algebraic_field(I)
    assert Domain.preprocess('QQ<I>') == QQ.algebraic_field(I)

    assert Domain.preprocess('Q<sqrt(2), I>') == QQ.algebraic_field(sqrt(2), I)
    assert Domain.preprocess(
        'QQ<sqrt(2), I>') == QQ.algebraic_field(sqrt(2), I)

    pytest.raises(OptionError, lambda: Domain.preprocess('abc'))

    assert Domain.preprocess('RR') == RR
    assert Domain.preprocess('RR_5') == RealField(prec=5)

    assert Domain.preprocess('CC') == CC
    assert Domain.preprocess('CC_5') == ComplexField(prec=5)

    pytest.raises(OptionError, lambda: Domain.preprocess(()))


def test_Domain_postprocess():
    pytest.raises(GeneratorsError, lambda: Domain.postprocess({'gens': (x, y),
                                                               'domain': ZZ.inject(y, z)}))

    pytest.raises(GeneratorsError, lambda: Domain.postprocess({'gens': (),
                                                               'domain': EX}))
    pytest.raises(GeneratorsError, lambda: Domain.postprocess({'domain': EX}))


def test_Split_preprocess():
    assert Split.preprocess(False) is False
    assert Split.preprocess(True) is True

    assert Split.preprocess(0) is False
    assert Split.preprocess(1) is True

    pytest.raises(OptionError, lambda: Split.preprocess(x))


def test_Split_postprocess():
    pytest.raises(NotImplementedError, lambda: Split.postprocess({'split': True}))


def test_Gaussian_preprocess():
    assert Gaussian.preprocess(False) is False
    assert Gaussian.preprocess(True) is True

    assert Gaussian.preprocess(0) is False
    assert Gaussian.preprocess(1) is True

    pytest.raises(OptionError, lambda: Gaussian.preprocess(x))


def test_Gaussian_postprocess():
    opt = {'gaussian': True}
    Gaussian.postprocess(opt)

    assert opt == {
        'gaussian': True,
        'extension': {I},
        'domain': QQ.algebraic_field(I),
    }


def test_Extension_preprocess():
    assert Extension.preprocess(True) is True
    assert Extension.preprocess(1) is True
    assert Extension.preprocess(False) is False

    assert Extension.preprocess([]) is None

    assert Extension.preprocess(sqrt(2)) == {sqrt(2)}
    assert Extension.preprocess([sqrt(2)]) == {sqrt(2)}

    assert Extension.preprocess([sqrt(2), I]) == {sqrt(2), I}


def test_Extension_postprocess():
    opt = {'extension': {sqrt(2)}}
    Extension.postprocess(opt)

    assert opt == {
        'extension': {sqrt(2)},
        'domain': QQ.algebraic_field(sqrt(2)),
    }

    opt = {'extension': True}
    Extension.postprocess(opt)

    assert opt == {'extension': True}


def test_Modulus_preprocess():
    assert Modulus.preprocess(23) == 23
    assert Modulus.preprocess(Integer(23)) == 23

    pytest.raises(OptionError, lambda: Modulus.preprocess(0))
    pytest.raises(OptionError, lambda: Modulus.preprocess(x))


def test_Modulus_postprocess():
    opt = {'modulus': 5}
    Modulus.postprocess(opt)

    assert opt == {
        'modulus': 5,
        'domain': FF(5),
    }

    opt = {'modulus': 5}
    Modulus.postprocess(opt)

    assert opt == {
        'modulus': 5,
        'domain': FF(5),
    }


def test_Strict_preprocess():
    assert Strict.preprocess(False) is False
    assert Strict.preprocess(True) is True

    assert Strict.preprocess(0) is False
    assert Strict.preprocess(1) is True

    pytest.raises(OptionError, lambda: Strict.preprocess(x))


def test_Strict_postprocess():
    opt = {'strict': True}
    Strict.postprocess(opt)

    assert opt == {'strict': True}


def test_Auto_preprocess():
    assert Auto.preprocess(False) is False
    assert Auto.preprocess(True) is True

    assert Auto.preprocess(0) is False
    assert Auto.preprocess(1) is True

    pytest.raises(OptionError, lambda: Auto.preprocess(x))


def test_Auto_postprocess():
    opt = {'auto': True}
    Auto.postprocess(opt)

    assert opt == {'auto': True}


def test_Frac_preprocess():
    assert Frac.preprocess(False) is False
    assert Frac.preprocess(True) is True

    assert Frac.preprocess(0) is False
    assert Frac.preprocess(1) is True

    pytest.raises(OptionError, lambda: Frac.preprocess(x))


def test_Frac_postprocess():
    opt = {'frac': True}
    Frac.postprocess(opt)

    assert opt == {'frac': True}


def test_Formal_preprocess():
    assert Formal.preprocess(False) is False
    assert Formal.preprocess(True) is True

    assert Formal.preprocess(0) is False
    assert Formal.preprocess(1) is True

    pytest.raises(OptionError, lambda: Formal.preprocess(x))


def test_Formal_postprocess():
    opt = {'formal': True}
    Formal.postprocess(opt)

    assert opt == {'formal': True}


def test_Polys_preprocess():
    assert Polys.preprocess(False) is False
    assert Polys.preprocess(True) is True

    assert Polys.preprocess(0) is False
    assert Polys.preprocess(1) is True

    pytest.raises(OptionError, lambda: Polys.preprocess(x))


def test_Polys_postprocess():
    opt = {'polys': True}
    Polys.postprocess(opt)

    assert opt == {'polys': True}


def test_Include_preprocess():
    assert Include.preprocess(False) is False
    assert Include.preprocess(True) is True

    assert Include.preprocess(0) is False
    assert Include.preprocess(1) is True

    pytest.raises(OptionError, lambda: Include.preprocess(x))


def test_Include_postprocess():
    opt = {'include': True}
    Include.postprocess(opt)

    assert opt == {'include': True}


def test_All_preprocess():
    assert All.preprocess(False) is False
    assert All.preprocess(True) is True

    assert All.preprocess(0) is False
    assert All.preprocess(1) is True

    pytest.raises(OptionError, lambda: All.preprocess(x))


def test_All_postprocess():
    opt = {'all': True}
    All.postprocess(opt)

    assert opt == {'all': True}


def test_Gen_preprocess():
    opt = {'gen': 'spam'}
    pytest.raises(OptionError, lambda: Gen.preprocess(opt))


def test_Gen_postprocess():
    opt = {'gen': x}
    Gen.postprocess(opt)

    assert opt == {'gen': x}


def test_Symbols_preprocess():
    pytest.raises(OptionError, lambda: Symbols.preprocess(x))


def test_Symbols_postprocess():
    opt = {'symbols': [x, y, z]}
    Symbols.postprocess(opt)

    assert opt == {'symbols': [x, y, z]}


def test_Method_preprocess():
    pytest.raises(OptionError, lambda: Method.preprocess(10))


def test_Method_postprocess():
    opt = {'method': 'f5b'}
    Method.postprocess(opt)

    assert opt == {'method': 'f5b'}


def test_allowed_flags():
    pytest.raises(OptionError, lambda: allowed_flags({'spam': True}, []))


def test_set_defaults():
    assert set_defaults({'defaults': None}) == {'defaults': None}
