"""
Utilities for our testing framework.
"""

from __future__ import print_function, division

import os
import sys
import subprocess
import re
from sympy.core.compatibility import PY3


IS_WINDOWS = (os.name == 'nt')


class Reporter(object):
    """
    Parent class for all reporters.
    """
    pass


class PyTestReporter(Reporter):
    """
    Py.test like reporter. Should produce output identical to py.test.
    """

    def __init__(self, verbose=False, tb="short", colors=True,
                 force_colors=False, split=None):
        self._verbose = verbose
        self._tb_style = tb
        self._colors = colors
        self._force_colors = force_colors
        self._xfailed = 0
        self._xpassed = []
        self._failed = []
        self._failed_doctest = []
        self._passed = 0
        self._skipped = 0
        self._exceptions = []
        self._terminal_width = None
        self._default_width = 80
        self._split = split

        # this tracks the x-position of the cursor (useful for positioning
        # things on the screen), without the need for any readline library:
        self._write_pos = 0
        self._line_wrap = False

    def root_dir(self, dir):
        self._root_dir = dir

    @property
    def terminal_width(self):
        if self._terminal_width is not None:
            return self._terminal_width

        def findout_terminal_width():
            if sys.platform == "win32":
                # Windows support is based on:
                #
                #  http://code.activestate.com/recipes/
                #  440694-determine-size-of-console-window-on-windows/

                from ctypes import windll, create_string_buffer

                h = windll.kernel32.GetStdHandle(-12)
                csbi = create_string_buffer(22)
                res = windll.kernel32.GetConsoleScreenBufferInfo(h, csbi)

                if res:
                    import struct
                    (_, _, _, _, _, left, _, right, _, _, _) = \
                        struct.unpack("hhhhHhhhhhh", csbi.raw)
                    return right - left
                else:
                    return self._default_width

            if hasattr(sys.stdout, 'isatty') and not sys.stdout.isatty():
                return self._default_width  # leave PIPEs alone

            try:
                process = subprocess.Popen(['stty', '-a'],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
                stdout = process.stdout.read()
                if PY3:
                    stdout = stdout.decode("utf-8")
            except (OSError, IOError):
                pass
            else:
                # We support the following output formats from stty:
                #
                # 1) Linux   -> columns 80
                # 2) OS X    -> 80 columns
                # 3) Solaris -> columns = 80

                re_linux = r"columns\s+(?P<columns>\d+);"
                re_osx = r"(?P<columns>\d+)\s*columns;"
                re_solaris = r"columns\s+=\s+(?P<columns>\d+);"

                for regex in (re_linux, re_osx, re_solaris):
                    match = re.search(regex, stdout)

                    if match is not None:
                        columns = match.group('columns')

                        try:
                            width = int(columns)
                        except ValueError:
                            pass
                        if width != 0:
                            return width

            return self._default_width

        width = findout_terminal_width()
        self._terminal_width = width

        return width

    def write(self, text, color="", align="left", width=None,
              force_colors=False):
        """
        Prints a text on the screen.

        It uses sys.stdout.write(), so no readline library is necessary.

        Parameters
        ==========

        color : choose from the colors below, "" means default color
        align : "left"/"right", "left" is a normal print, "right" is aligned on
                the right-hand side of the screen, filled with spaces if
                necessary
        width : the screen width

        """
        color_templates = (
            ("Black", "0;30"),
            ("Red", "0;31"),
            ("Green", "0;32"),
            ("Brown", "0;33"),
            ("Blue", "0;34"),
            ("Purple", "0;35"),
            ("Cyan", "0;36"),
            ("LightGray", "0;37"),
            ("DarkGray", "1;30"),
            ("LightRed", "1;31"),
            ("LightGreen", "1;32"),
            ("Yellow", "1;33"),
            ("LightBlue", "1;34"),
            ("LightPurple", "1;35"),
            ("LightCyan", "1;36"),
            ("White", "1;37"),
        )

        colors = {}

        for name, value in color_templates:
            colors[name] = value
        c_normal = '\033[0m'
        c_color = '\033[%sm'

        if width is None:
            width = self.terminal_width

        if align == "right":
            if self._write_pos + len(text) > width:
                # we don't fit on the current line, create a new line
                self.write("\n")
            self.write(" "*(width - self._write_pos - len(text)))

        if not self._force_colors and hasattr(sys.stdout, 'isatty') and not \
                sys.stdout.isatty():
            # the stdout is not a terminal, this for example happens if the
            # output is piped to less, e.g. "bin/test | less". In this case,
            # the terminal control sequences would be printed verbatim, so
            # don't use any colors.
            color = ""
        elif sys.platform == "win32":
            # Windows consoles don't support ANSI escape sequences
            color = ""
        elif not self._colors:
            color = ""

        if self._line_wrap:
            if text[0] != "\n":
                sys.stdout.write("\n")

        # Avoid UnicodeEncodeError when printing out test failures
        if PY3 and IS_WINDOWS:
            text = text.encode('raw_unicode_escape').decode('utf8', 'ignore')
        elif PY3 and not sys.stdout.encoding.lower().startswith('utf'):
            text = text.encode(sys.stdout.encoding, 'backslashreplace'
                              ).decode(sys.stdout.encoding)

        if color == "":
            sys.stdout.write(text)
        else:
            sys.stdout.write("%s%s%s" %
                (c_color % colors[color], text, c_normal))
        sys.stdout.flush()
        l = text.rfind("\n")
        if l == -1:
            self._write_pos += len(text)
        else:
            self._write_pos = len(text) - l - 1
        self._line_wrap = self._write_pos >= width
        self._write_pos %= width

    def write_center(self, text, delim="="):
        width = self.terminal_width
        if text != "":
            text = " %s " % text
        idx = (width - len(text)) // 2
        t = delim*idx + text + delim*(width - idx - len(text))
        self.write(t + "\n")

    def write_exception(self, e, val, tb):
        t = traceback.extract_tb(tb)
        # remove the first item, as that is always runtests.py
        t = t[1:]
        t = traceback.format_list(t)
        self.write("".join(t))
        t = traceback.format_exception_only(e, val)
        self.write("".join(t))

    def start(self, seed=None, msg="test process starts"):
        self.write_center(msg)
        executable = sys.executable
        v = tuple(sys.version_info)
        python_version = "%s.%s.%s-%s-%s" % v
        implementation = platform.python_implementation()
        if implementation == 'PyPy':
            implementation += " %s.%s.%s-%s-%s" % sys.pypy_version_info
        self.write("executable:         %s  (%s) [%s]\n" %
            (executable, python_version, implementation))
        from .misc import ARCH
        self.write("architecture:       %s\n" % ARCH)
        from sympy.core.cache import USE_CACHE
        self.write("cache:              %s\n" % USE_CACHE)
        from sympy.core.compatibility import GROUND_TYPES, HAS_GMPY
        version = ''
        if GROUND_TYPES =='gmpy':
            if HAS_GMPY == 1:
                import gmpy
            elif HAS_GMPY == 2:
                import gmpy2 as gmpy
            version = gmpy.version()
        self.write("ground types:       %s %s\n" % (GROUND_TYPES, version))
        if seed is not None:
            self.write("random seed:        %d\n" % seed)
        from .misc import HASH_RANDOMIZATION
        self.write("hash randomization: ")
        hash_seed = os.getenv("PYTHONHASHSEED") or '0'
        if HASH_RANDOMIZATION and (hash_seed == "random" or int(hash_seed)):
            self.write("on (PYTHONHASHSEED=%s)\n" % hash_seed)
        else:
            self.write("off\n")
        if self._split:
            self.write("split:              %s\n" % self._split)
        self.write('\n')
        self._t_start = clock()

    def finish(self):
        self._t_end = clock()
        self.write("\n")
        global text, linelen
        text = "tests finished: %d passed, " % self._passed
        linelen = len(text)

        def add_text(mytext):
            global text, linelen
            """Break new text if too long."""
            if linelen + len(mytext) > self.terminal_width:
                text += '\n'
                linelen = 0
            text += mytext
            linelen += len(mytext)

        if len(self._failed) > 0:
            add_text("%d failed, " % len(self._failed))
        if len(self._failed_doctest) > 0:
            add_text("%d failed, " % len(self._failed_doctest))
        if self._skipped > 0:
            add_text("%d skipped, " % self._skipped)
        if self._xfailed > 0:
            add_text("%d expected to fail, " % self._xfailed)
        if len(self._xpassed) > 0:
            add_text("%d expected to fail but passed, " % len(self._xpassed))
        if len(self._exceptions) > 0:
            add_text("%d exceptions, " % len(self._exceptions))
        add_text("in %.2f seconds" % (self._t_end - self._t_start))

        if len(self._xpassed) > 0:
            self.write_center("xpassed tests", "_")
            for e in self._xpassed:
                self.write("%s: %s\n" % (e[0], e[1]))
            self.write("\n")

        if self._tb_style != "no" and len(self._exceptions) > 0:
            for e in self._exceptions:
                filename, f, (t, val, tb) = e
                self.write_center("", "_")
                if f is None:
                    s = "%s" % filename
                else:
                    s = "%s:%s" % (filename, f.__name__)
                self.write_center(s, "_")
                self.write_exception(t, val, tb)
            self.write("\n")

        if self._tb_style != "no" and len(self._failed) > 0:
            for e in self._failed:
                filename, f, (t, val, tb) = e
                self.write_center("", "_")
                self.write_center("%s:%s" % (filename, f.__name__), "_")
                self.write_exception(t, val, tb)
            self.write("\n")

        if self._tb_style != "no" and len(self._failed_doctest) > 0:
            for e in self._failed_doctest:
                filename, msg = e
                self.write_center("", "_")
                self.write_center("%s" % filename, "_")
                self.write(msg)
            self.write("\n")

        self.write_center(text)
        ok = len(self._failed) == 0 and len(self._exceptions) == 0 and \
            len(self._failed_doctest) == 0
        if not ok:
            self.write("DO *NOT* COMMIT!\n")
        return ok

    def entering_filename(self, filename, n):
        rel_name = filename[len(self._root_dir) + 1:]
        self._active_file = rel_name
        self._active_file_error = False
        self.write(rel_name)
        self.write("[%d] " % n)

    def leaving_filename(self):
        self.write(" ")
        if self._active_file_error:
            self.write("[FAIL]", "Red", align="right")
        else:
            self.write("[OK]", "Green", align="right")
        self.write("\n")
        if self._verbose:
            self.write("\n")

    def entering_test(self, f):
        self._active_f = f
        if self._verbose:
            self.write("\n" + f.__name__ + " ")

    def test_xfail(self):
        self._xfailed += 1
        self.write("f", "Green")

    def test_xpass(self, v):
        message = str(v)
        self._xpassed.append((self._active_file, message))
        self.write("X", "Green")

    def test_fail(self, exc_info):
        self._failed.append((self._active_file, self._active_f, exc_info))
        self.write("F", "Red")
        self._active_file_error = True

    def doctest_fail(self, name, error_msg):
        # the first line contains "******", remove it:
        error_msg = "\n".join(error_msg.split("\n")[1:])
        self._failed_doctest.append((name, error_msg))
        self.write("F", "Red")
        self._active_file_error = True

    def test_pass(self, char="."):
        self._passed += 1
        if self._verbose:
            self.write("ok", "Green")
        else:
            self.write(char, "Green")

    def test_skip(self, v=None):
        char = "s"
        self._skipped += 1
        if v is not None:
            message = str(v)
            if message == "KeyboardInterrupt":
                char = "K"
            elif message == "Timeout":
                char = "T"
            elif message == "Slow":
                char = "w"
        self.write(char, "Blue")
        if self._verbose:
            self.write(" - ", "Blue")
            if v is not None:
                self.write(message, "Blue")

    def test_exception(self, exc_info):
        self._exceptions.append((self._active_file, self._active_f, exc_info))
        self.write("E", "Red")
        self._active_file_error = True

    def import_error(self, filename, exc_info):
        self._exceptions.append((filename, None, exc_info))
        rel_name = filename[len(self._root_dir) + 1:]
        self.write(rel_name)
        self.write("[?]   Failed to import", "Red")
        self.write(" ")
        self.write("[FAIL]", "Red", align="right")
        self.write("\n")
