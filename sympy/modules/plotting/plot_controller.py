from pyglet.window import key
from pyglet.window.mouse import LEFT

class PlotController(object):
    
    normal_mouse_sensitivity = 4.0
    modified_mouse_sensitivity = 1.0

    normal_key_sensitivity = 16.0
    modified_key_sensitivity = 4.0

    keymap = {
                key.LEFT:'left',
                key.A:'left',
                
                key.RIGHT:'right',
                key.D:'right',
                
                key.UP:'up',
                key.W:'up',
                
                key.DOWN:'down',
                key.S:'down',
                
                key.NUM_ADD:'zoom_in',
                key.PAGEUP:'zoom_in',
                key.R:'zoom_in',

                key.NUM_SUBTRACT:'zoom_out',
                key.PAGEDOWN:'zoom_out',
                key.F:'zoom_out',
                
                key.RSHIFT:'modify_sensitivity',
                key.LSHIFT:'modify_sensitivity',
             }

    def __init__(self, window):
        self.action = {

                'left':False,
                'right':False,
                'up':False,
                'down':False,
                'zoom_in':False, 
                'zoom_out':False,
                'modify_sensitivity':False,

            }
        self.window = window
        
    def update(self, dt):
        z = 0
        if self.action['zoom_out']: z -= 1
        if self.action['zoom_in']: z += 1
        if z != 0:
            self.window.camera.zoom_relative(z/10.0, self.get_key_sensitivity())        
        
        dx, dy = 0, 0
        if self.action['left']: dx -= 1
        if self.action['right']: dx += 1
        if self.action['down']: dy -= 1
        if self.action['up']: dy += 1

        if dx != 0 or dy != 0:
            dx = float(dx) * dt * 100.0 # yep, a magic number
            dy = float(dy) * dt * 100.0

            p1 = (self.window.width/2, self.window.height/2)
            p2 = ( p1[0] + dx, p1[1] + dy )
            self.window.camera.spherical_rotate(p1, p2, self.get_key_sensitivity())

        return True

    def get_mouse_sensitivity(self):
        if self.action['modify_sensitivity']:
            return self.modified_mouse_sensitivity
        else:
            return self.normal_mouse_sensitivity

    def get_key_sensitivity(self):
        if self.action['modify_sensitivity']:
            return self.modified_key_sensitivity
        else:
            return self.normal_key_sensitivity

    def on_key_press(self, symbol, modifiers):
        if symbol in self.keymap:
            self.action[self.keymap[symbol]] = True

    def on_key_release(self, symbol, modifiers):
        if symbol in self.keymap:
            self.action[self.keymap[symbol]] = False

    def on_mouse_drag(self, x, y, dx, dy, buttons, modifiers):
        if buttons & LEFT:
            self.window.camera.spherical_rotate((x-dx,y-dy),(x,y), self.get_mouse_sensitivity())

    def on_mouse_scroll(self, x, y, dx, dy):
        self.window.camera.zoom_relative(dy, self.get_mouse_sensitivity())
