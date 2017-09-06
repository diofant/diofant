"""ASCII-ART 2D pretty-printer"""

from .pretty import (pretty, pretty_print, pprint,  # noqa: F401
                     pprint_use_unicode, pprint_try_use_unicode)

# if unicode output is available -- let's use it
pprint_try_use_unicode()
