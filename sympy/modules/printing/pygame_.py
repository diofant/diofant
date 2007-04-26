#install TeX and these Debian packages: python-pygame, python-pexpect, dvipng
#

import tempfile
from sympy.modules.printing import print_latex

def print_pygame(st):

    try:
        import pygame
    except ImportError:
        print "Pygame is not installed. In Debian, install the " \
            "python-pygame package."
        return

    st = print_latex(st)

    from pygame import QUIT, KEYDOWN, K_ESCAPE, K_q
    pygame.font.init()
    size = 640, 240
    screen = pygame.display.set_mode(size)
    screen.fill((255, 255, 255))
    font = pygame.font.Font(None, 24)
    text = font.render(st, True, (0, 0, 0))
    textpos = text.get_rect(centerx=screen.get_width()/2)
    screen.blit(text, textpos)
    pygame.display.flip()

    image = tex2png(st,pygame)
    imagepos = image.get_rect(centerx=screen.get_width()/2).move((0,30))
    screen.blit(image, imagepos)
    pygame.display.flip()

    while 1:
        for event in pygame.event.get():
            if event.type == QUIT:
                return
            elif event.type == KEYDOWN and event.key == K_ESCAPE:
                return
            elif event.type == KEYDOWN and event.key == K_q:
                return

def tex2png(eq,pygame):
    #http://www.fauskes.net/nb/htmleqII/
    import os, sys

    import pexpect

    # Include your favourite plain TeX packages and commands here
    tex_preamble = r"""\nopagenumbers
"""

    x = tempfile.mktemp()
    tmp1 = '%s.tex'%x

    # create a LaTeX document and insert equations
    f = open(tmp1,'w')
    f.write(tex_preamble)
    f.write(r"""$%s$
\vfill
\eject
""" % eq)
    # end LaTeX document    
    f.write(r'\bye')
    f.close()

    # compile LaTeX document. A DVI file is created
    cwd = os.getcwd()
    os.chdir("/tmp")
    pexpect.run('tex %s' % tmp1)

    # Run dvipng on the generated DVI file. Use tight bounding box. 
    # Magnification is set to 1200
    # currently, the dvipng is broken on debian.....
    cmd = "dvipng -T tight -x 1728 -z 9 -bg transparent " \
    + "-o %s.png %s.dvi" % (x,x)
    pexpect.run(cmd) 
    image = pygame.image.load("%s.png" % x)

    #remove temporary files
    os.remove("%s.tex" % x)
    os.remove("%s.dvi" % x)
    os.remove("%s.log" % x)
    os.remove("%s.png" % x)
    os.chdir(cwd)

    return image
