import sys
sys.path.append(".")

from sym import hashing

def test_basics():
    h1=hashing.mhash()
    h1.addint(-1)
    h1.addint(-2)

    h2=hashing.mhash()
    h2.addint(-1)
    h2.addint(-2)
    assert h1.value == h2.value

    h3=hashing.mhash()
    h3.addint(-2)
    h3.addint(-1)
    assert h1.value != h3.value
