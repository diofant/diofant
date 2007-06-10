
from math import ceil

class BitField(object):
    """Implementation of bit-field data structure based on Python lists
       of 32-bit integers, although its capapcity is not aligned to this
       boundary. This class exports container facilities so it can be
       used in iterable contexts but for simplicity negative indexes
       are not allowed.

       >>> bits = BitField(123)
       >>> bits[0] = True
       >>> bits[17] = True

       >>> bits[3]
       False
       >>> bits[17]
       True
       >>> len(bits)
       123

       >>> len([ bit for bit in bits if bit ])
       2

    """

    def __init__(self, size, default=False):
        self.size, self.cells = size, int(ceil(size/32.0))

        if default == True:
            self.bits = [0xFFFFFFFF] * self.cells
        else:
            self.bits = [0] * self.cells

    def __len__(self):
        return self.size

    def __getitem__(self, key):
        if 0 <= key and key < self.size:
            return bool(self.bits[key>>5] & (1<<(key-((key>>5)<<5))))
        else:
            raise ValueError("Index out of bounds")

    def __setitem__(self, key, value):
        if 0 <= key and key < self.size:
            cell = key >> 5
            mask = 1 << (key-(cell<<5))

            if value == True:
                self.bits[cell] |= mask
            else:
                self.bits[cell] &= 0xFFFFFFFF ^ mask
        else:
            raise ValueError("Index out of bounds")

    def __iter__(self):
        for i in range(self.cells-1):
            for j in range(32):
                yield bool(self.bits[i] & (1<<j))

        for i in range(self.size - ((self.size>>5)<<5)):
            yield bool(self.bits[-1] & (1<<i))

