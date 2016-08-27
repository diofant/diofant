#!/usr/bin/env python

DESCRIPTION = """
Runs all the examples for testing purposes and reports successes and failures
to stderr.  An example is marked successful if the running thread does not
throw an exception, for threaded examples, such as plotting, one needs to
check the stderr messages as well.
"""

EPILOG = """
Example Usage:
   When no examples fail:
     $ ./all.py > out
     SUCCESSFUL:
       - beginner.basic
       [...]
     NO FAILED EXAMPLES
     $

   When examples fail:
     $ ./all.py -w > out
     Traceback (most recent call last):
       File "./all.py", line 111, in run_examples
     [...]
     SUCCESSFUL:
       - beginner.basic
       [...]
     FAILED:
       - intermediate.mplot2D
       [...]
     $

   Obviously, we want to achieve the first result.
"""

import imp
import optparse
import os
import sys
import traceback

# add local diofant to the module path
this_file = os.path.abspath(__file__)
diofant_dir = os.path.join(os.path.dirname(this_file), "..")
diofant_dir = os.path.normpath(diofant_dir)
sys.path.insert(0, diofant_dir)

TERMINAL_EXAMPLES = [
    "intermediate.partial_differential_eqs",
    "intermediate.trees",
    "intermediate.vandermonde",
    "advanced.curvilinear_coordinates",
    "advanced.fem",
    "advanced.gibbs_phenomenon",
    "advanced.relativity",
]

WINDOWED_EXAMPLES = [
    "advanced.autowrap_ufuncify",
]

EXAMPLE_DIR = os.path.dirname(__file__)


def __import__(name, globals=None, locals=None, fromlist=None):
    """An alternative to the import function so that we can import
    modules defined as strings.

    This code was taken from: http://docs.python.org/lib/examples-imp.html
    """
    # Fast path: see if the module has already been imported.
    try:
        return sys.modules[name]
    except KeyError:
        pass

    # If any of the following calls raises an exception,
    # there's a problem we can't handle -- let the caller handle it.
    module_name = name.split('.')[-1]
    module_path = os.path.join(EXAMPLE_DIR, *name.split('.')[:-1])

    fp, pathname, description = imp.find_module(module_name, [module_path])

    try:
        return imp.load_module(module_name, fp, pathname, description)
    finally:
        # Since we may exit via an exception, close fp explicitly.
        if fp:
            fp.close()


def load_example_module(example):
    """Loads modules based upon the given package name"""
    mod = __import__(example)
    return mod


def run_examples(windowed=False, quiet=False, summary=True):
    """Run all examples in the list of modules.

    Returns a boolean value indicating whether all the examples were
    successful.
    """
    successes = []
    failures = []
    examples = TERMINAL_EXAMPLES
    if windowed:
        examples += WINDOWED_EXAMPLES

    for example in examples:
        if run_example(example, quiet):
            successes.append(example)
        else:
            failures.append(example)

    if summary:
        show_summary(successes, failures, quiet)

    return len(failures) == 0


def run_example(example, quiet=False):
    """Run a specific example.

    Returns a boolean value indicating whether the example was successful.
    """
    if quiet:
        print(example + " " * (72 - len(example)), end='')
    else:
        print("=" * 79)
        print("Running: ", example)

    try:
        mod = load_example_module(example)
        if quiet:
            suppress_output(mod.main)
            print("[PASS]")
        else:
            mod.main()
        return True
    except:
        if quiet:
            print("[FAIL]")
        traceback.print_exc()
        return False


class DummyFile(object):
    def write(self, x):
        pass


def suppress_output(fn):
    """Suppresses the output of fn on sys.stdout."""
    save_stdout = sys.stdout
    try:
        sys.stdout = DummyFile()
        fn()
    finally:
        sys.stdout = save_stdout


def show_summary(successes, failures, quiet=False):
    """Shows a summary detailing which examples were successful and which failed."""
    if quiet:
        print("-" * 79)
        if failures:
            print("FAILED:")
            for example in failures:
                print("  " + example)
        else:
            print("ALL EXAMPLES PASSED")
    else:
        if successes:
            print("SUCCESSFUL: ", file=sys.stderr)
            for example in successes:
                print("  -", example, file=sys.stderr)
        else:
            print("NO SUCCESSFUL EXAMPLES", file=sys.stderr)

        if failures:
            print("FAILED: ", file=sys.stderr)
            for example in failures:
                print("  -", example, file=sys.stderr)
        else:
            print("NO FAILED EXAMPLES", file=sys.stderr)


def main(*args, **kws):
    """Main script runner"""
    parser = optparse.OptionParser()
    parser.add_option('-w', '--windowed', action="store_true", dest="windowed",
        help="also run examples requiring windowed environment")
    parser.add_option('-q', '--quiet', action="store_true", dest="quiet",
        help="runs examples in 'quiet mode' suppressing example output and \
              showing simple status messages.")
    parser.add_option('--no-summary', action="store_true", dest="no_summary",
        help="hides the summary at the end of testing the examples")

    (options, _) = parser.parse_args()

    return 0 if run_examples(windowed=options.windowed, quiet=options.quiet,
                             summary=not options.no_summary) else 1


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
