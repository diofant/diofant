def print_pygame(st):
    try:
        import pygame
    except ImportError:
        print "Pygame is not installed. In Debian, install the " \
            "python-pygame package."
        return
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

    while 1:
        for event in pygame.event.get():
            if event.type == QUIT:
                return
            elif event.type == KEYDOWN and event.key == K_ESCAPE:
                return
            elif event.type == KEYDOWN and event.key == K_q:
                return

        #screen.blit(text, textpos)
        #pygame.display.flip()
