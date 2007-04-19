#!/usr/bin/env python
"""Distutils based setup script for Sympy.

This uses Distutils (http://python.org/sigs/distutils-sig/) the standard
python mechanism for installing packages. For the easiest installation
just type the command (you'll probably need root privileges for that):

    python setup.py install

This will install the library in the default location. For instructions on
how to customize the install procedure read the output of:

    python setup.py --help install

In addition, there are some other commands:

    python setup.py bdist_dpkg -> will make a deb package in the parent diretory
    python setup.py clean -> will clean all trash (*.pyc and stuff)
    python setup.py test  -> will run the complete test suite
    python setup.py test_core -> will run only tests concerning core features
    python setup.py test_doc -> will run tests on the examples of the documentation
    
To get a full list of avaiable commands, read the output of: 

    python setup.py --help-commands
    
Or, if all else fails, feel free to write to the sympy list at
sympy@googlegroups.com and ask for help.
"""

from distutils.core import setup
from distutils.core import Command
import sys
import os

import sympy

# Make sure I have the right Python version.
if sys.version_info[1] < 4:
    print "Sympy requires Python 2.4 or newer.  Python %d.%d detected" % \
          sys.version_info[:2]
    sys.exit(-1)

class test_sympy(Command):
    """Runs all tests under the tests/ folder
    """
    
    description = "Automatically run the test suite for Sympy."
    user_options = []  # distutils complains if this is not here.


    def initialize_options(self):  # distutils wants this
        pass
    
    def finalize_options(self):    # this too
        pass
    
    def run(self):
        try:
            import py
        except ImportError:
            print """In order to run the tests, you need codespeak's py.lib
            web page: http://codespeak.net/py/dist/
            If you are on debian systems, the package is named python-codespeak-lib
            """
            sys.exit(-1)
        py.test.cmdline.main(args=["tests"])
        test_sympy_doc.run_doctest() # run also the doc test suite
        
class bdist_dpkg(Command):
    """Make a nice .deb package
    """
    
    description = "Make a deb package using dpkg (debuild)"
    user_options = []  # distutils complains if this is not here.

    def initialize_options(self):  # distutils wants this
        pass
    
    def finalize_options(self):    # this too
        pass

    def run(self):
        """
        debian/changelog contains a version like this:

        0.4~pre+svn739-1

        This method parses it, then checkouts the svn revision as directed (739
        in this example), but applies the current top svn debian dir to it, and
        executes "debuild" in that temporary directory.
        """
        import os
        def get_changelog_version_revision():
            """Reads the first line in changelog, parses 0.4~pre+svn739-1 and
            returns ("0.4~pre",739,1)
            """
            l = file("debian/changelog").readline()
            import re
            m = re.match("sympy \((\S+)\+svn(\d+)\-(\d+)\) ",l)
            if m:
                g = m.groups()
                if len(g) == 3:
                    #version, svn revision, debian revision
                    #('0.4~pre', '739', '1') 
                    v, r, dr = g
                    return v, int(r), int(dr)
            print l
            raise "Don't understant the syntax in changelog"
        version,revision,drevision = get_changelog_version_revision()
        os.system("mkdir -p dist")
        tmpdir = "sympy-%s+svn%d" % (version, revision)
        print "exporting svn (%d) to dist/%s" % (revision,tmpdir)
        os.system("svn -q export -r %d " % revision +
                "http://sympy.googlecode.com/svn/trunk/ dist/%s" % tmpdir)  
        os.system("rm -rf dist/%s/debian" % tmpdir)
        print "creating dist/sympy_%s+svn%d.orig.tar.gz" % (version, revision)
        os.system("cd dist; tar zcf sympy_%s+svn%d.orig.tar.gz %s" \
                %(version, revision, tmpdir))
        print "creating the deb package"
        os.system("cp -a debian dist/%s/debian" % tmpdir)
        os.system("rm -rf dist/%s/debian/.svn" % tmpdir)
        #os.system("cd dist/%s; debuild -sa -us -uc" % tmpdir)
        os.system("cd dist/%s; debuild" % tmpdir)
        os.system("rm -rf dist/%s" % tmpdir)
        print "-"*50
        print "Done. Files genereated in the dist/ directory"
    
    def run_old(self):
        """
        Copies the current local svn copy to the dist/sympy-svn739,
        creates dist/sympy_0.4-pre+svn739.orig.tar.gz from that and then
        a debian package.

        debian/changelog contains a version like this:

        0.4-pre+svn739-1

        where the last "-1" is a debian version and this is the only place that
        contains this number. The "0.4-pre" should be the same as
        sympy.__init__ and svn739 should be the same as the current svn number.

        This method checks that 0.4-pre+svn739 in changelog is consistent,
        otherwise refuses to continue.
        """
        import os
        def get_revision():
            """Returns the current svn revision number as an int."""
            fin, fout = os.popen2("svn --non-interactive info | grep Revision | cut -c 11-")
            rev = fout.readlines()
            return int(rev[0])
        def get_changelog_version_revision():
            """Reads the first line in changelog, parses 0.4-pre+svn739-1 and
            returns ("0.4-pre",739,1)
            """
            l = file("debian/changelog").readline()
            import re
            m = re.match("sympy \((\S+)\+svn(\d+)\-(\d+)\) ",l)
            if m:
                g = m.groups()
                if len(g) == 3:
                    #version, svn revision, debian revision
                    #('0.4-pre', '739', '1') 
                    v, r, dr = g
                    return v, int(r), int(dr)
            print l
            raise "Don't understant the syntax in changelog"
        v,r,dr = get_changelog_version_revision()
        revision=get_revision()
        version = sympy.__version__
        pos = version.find("-")
        if pos != -1:
            #change "-" to "~" in the Debian package name run this to see why:
#$ dpkg --compare-versions 0.4-pre+svn756-1 lt 0.4 && echo true || echo false
#$ dpkg --compare-versions 0.4~pre+svn756-1 lt 0.4 && echo true || echo false
#$ dpkg --compare-versions 0.3+svn756-1 lt 0.4 && echo true || echo false
            version = version[:pos]+"~"+version[pos+1:]
        if version != v or revision != r:
            print "The version/revision in debian/changelog (%s-svn%d) is " \
                "inconsistent\nwith the sympy.__version__ (%s) and svn " \
                "revision (%d)" % (v, r, sympy.__version__, revision)
            return
        os.system("mkdir -p dist")
        tmpdir = "sympy-%s+svn%d" % (version, revision)
        print "exporting svn to dist/%s" % tmpdir
        os.system("svn -q export . dist/%s" % tmpdir)  
        os.system("rm -rf dist/%s/debian" % tmpdir)
        print "creating dist/sympy_%s+svn%d.orig.tar.gz" \
                % (version, revision)
        os.system("cd dist; tar zcf sympy_%s+svn%d.orig.tar.gz %s" \
                %(version, revision, tmpdir))
        print "creating the deb package"
        os.system("cp -a debian dist/%s/debian" % tmpdir)
        os.system("rm -rf dist/%s/debian/.svn" % tmpdir)
        #os.system("cd dist/%s; debuild -sa -us -uc" % tmpdir)
        os.system("cd dist/%s; debuild" % tmpdir)
        os.system("rm -rf dist/%s" % tmpdir)
        print "-"*50
        print "Done. Files genereated in the dist/ directory"

class clean(Command):
    """Cleans *.pyc and debian trashs, so you should get the same copy as 
    is in the svn.
    """
    
    description = "Clean everything"
    user_options = [("all","a","the same")]  

    def initialize_options(self):  
        self.all = None
    
    def finalize_options(self):   
        pass

    def run(self):
        import os
        os.system("py.cleanup")
        os.system("rm -f python-build-stamp-2.4")
        os.system("rm -f MANIFEST")
        os.system("rm -rf build")
        os.system("rm -rf dist")

class test_sympy_core(Command):
    """Run only the tests concerning features of sympy.core.
    It's a lot faster than running the complete test suite.
    """
    
    description = "Automatically run the core test suite for Sympy."
    user_options = []  # distutils complains if this is not here.

    tests_to_run = ["tests/test_arit.py", "tests/test_basic.py", 
                   "tests/test_diff.py", "tests/test_equal.py", 
                   "tests/test_eval.py", "tests/test_evalf.py", 
                   "tests/test_functions.py", "tests/test_hashing.py", 
                   "tests/test_numbers.py", "tests/test_series.py", 
                   "tests/test_str.py", "tests/test_subs.py", 
                   "tests/test_symbol.py"
                   ]
    
    def initialize_options(self):  # distutils wants this
        pass
    
    def finalize_options(self):    # this too
        pass

        
    def run(self):
        try:
            import py
        except ImportError:
            print """In order to run the tests, you need codespeak's py.lib
            web page: http://codespeak.net/py/dist/
            If you are on debian systems, the package is named python-codespeak-lib
            """
            sys.exit(-1)
        py.test.cmdline.main(args=self.tests_to_run)

class test_sympy_dpkg(Command):
    
    description = "Run tests for the deb package"
    user_options = []  # distutils complains if this is not here.
    
    def initialize_options(self):  # distutils wants this
        pass
    
    def finalize_options(self):    # this too
        pass
    
    def run(self):
        import os
        from glob import glob
        g = glob("dist/*.changes")
        assert len(g) == 1
        changes = g[0]
        g = glob("dist/*.dsc")
        assert len(g) == 1
        dsc = g[0]
        g = glob("dist/*.deb")
        assert len(g) == 1
        deb = g[0]
        print "testing this package:"
        print "  ",dsc
        print "  ",changes
        print "  ",deb
        print
        print "running lintian & linda..."
        os.system("lintian -i %s" % changes)
        os.system("linda -i %s" % changes)
        print 'running pbuilder (please run "sudo pbuilder update" yourself)...'
        os.system("sudo pbuilder build %s" % dsc)
        print "running piuparts"
        os.system("sudo piuparts -p %s" % deb)
        print "Done, see the above output for any errors."
        
class test_sympy_doc(Command):
    
    description = "Run the tests for the examples in the documentation"
    user_options = []  # distutils complains if this is not here.
    
    def initialize_options(self):  # distutils wants this
        pass
    
    def finalize_options(self):    # this too
        pass
    
    @staticmethod
    def run_doctest():
        
        import unittest
        import doctest
        
        import glob

        print "Testing docstrings."

        files = glob.glob('sympy/*/*.py') + glob.glob('sympy/modules/*/*.py')
        #make it work on Windows too:
        files = [f.replace("\\","/") for f in files]
        
        # files without doctests or that don't work
        files.remove('sympy/modules/printing/pygame_.py')
        files.remove('sympy/modules/printing/pretty.py') # see issue 53

        
        #testing for optional libraries
        try:
            import libxslt
        except ImportError:
            #remove tests that make use of libxslt1
            files.remove('sympy/modules/printing/latex.py')
            files.remove('sympy/modules/printing/__init__.py')
        try:
            import matplotlib
        except ImportError:
            files.remove('sympy/modules/graphing.py')

        modules = []
        for x in files:
            if len(x) > 12 and x[-11:] == '__init__.py':
                x = x.replace('/__init__', '') 
                print x
            modules.append(x.replace('/', '.')[:-3])
            #put . as separator and strip the extension (.py)

        modules.append('sympy')
        
        suite = unittest.TestSuite()
        for mod in modules:
            suite.addTest(doctest.DocTestSuite(mod))
            
        runner = unittest.TextTestRunner()
        runner.run(suite)

    def run(self):
        self.run_doctest()
        
import sympy

setup(
      name = 'sympy', 
      version = sympy.__version__, 
      description = 'Computer algebra system (CAS) in Python', 
      url = 'http://code.google.com/p/sympy', 
      packages = ['sympy', 'sympy.core', 'sympy.modules',
                  'sympy.modules.mathml', 'sympy.modules.printing'],
      package_data = {'sympy.modules.mathml' : ['data/mmlctop.xsl', 
                                                'data/mmltex.xsl',
                                                'data/simple_mmlctop.xsl' ]},
      scripts = ['bin/isympy'],
      ext_modules = [],
      data_files = [('share/man/man1', ['doc/man/isympy.1'])],
      cmdclass    = {'test': test_sympy, 
                     'test_core' : test_sympy_core,
                     'test_doc' : test_sympy_doc,
                     'bdist_dpkg' : bdist_dpkg, 
                     'test_dpkg' : test_sympy_dpkg,
                     'clean' : clean, 
                     },
      )

