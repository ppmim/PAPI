"""
   Main routines for the loggin function.
"""

import logging


"""logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s %(levelname)-8s %(module)s:%(lineno)d: %(message)s',
    datefmt='%a, %d %b %Y %H:%M:%S',
    filename='/tmp/panicTool.log',
    filemode='w')
#console = logging.StreamHandler()
#console.setLevel(logging.DEBUG)
#formatter = logging.Formatter('%(filename)s %(name)-12s: %(levelname)-8s %(message)s')
# tell the handler to use this format
#console.setFormatter(formatter)
#logging.getLogger('').addHandler(console)
log = logging.getLogger('panicQL.log')
log.setLevel(logging.INFO)
"""

class ColorFormatter(logging.Formatter):
 
    def color(self, level=None):
        codes = {\
            None:       (0,   0),
            'INFO':     (0,  33), # yellow
            'DEBUG':    (1,  32), # green
            'WARNING':  (1,  34), # blue
            'ERROR':    (1,  31), # red
            'CRITICAL': (1, 101), # negro, fondo rojo
            }
        return (chr(27)+'[%d;%dm') % codes[level]
 
    def format(self, record):
        retval = logging.Formatter.format(self, record)
        return self.color(record.levelname) + retval + self.color()

# Console
console = logging.StreamHandler()
console.setLevel(logging.INFO)
console.setFormatter(ColorFormatter('    [%(name)s]: %(asctime)s %(levelname)-8s %(module)s:%(lineno)d: %(message)s'))
logging.getLogger('panicQL').addHandler(console)
# File
file_hd = logging.FileHandler("/tmp/panicQL.log")
console.setLevel(logging.DEBUG)
formatter = logging.Formatter('[%(name)s]: %(asctime)s %(levelname)-8s %(module)s:%(lineno)d: %(message)s')
file_hd.setFormatter(formatter)
logging.getLogger('panicQL').addHandler(file_hd)

logging.getLogger('panicQL').setLevel(logging.DEBUG)
log = logging.getLogger('panicQL')


## Test
"""
log.debug('mensaje para depuracion')
log.info('informacion')
log.warning('el que avisa no es traidor')
log.error('un errorcillo')
log.critical('y la liaste parda')
"""
"""log = logging.getLogger('panicQL')
hdlr = logging.FileHandler('/tmp/panicQL.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
log.addHandler(hdlr)
log.setLevel(logging.INFO)
log.debug("Que pasa !!!!")
"""



def initLog(log_filename, log_level):
    return
    logging.basicConfig(level=log_level,
                        format='%(asctime)s %(levelname)-8s %(module)s:%(lineno)d: %(message)s',
                        datefmt='%a, %d %b %Y %H:%M:%S',
                        filename=log_filename,
                        filemode='w')
    




