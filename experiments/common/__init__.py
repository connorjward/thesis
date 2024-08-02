from . import petsclog
from . import roofline

def timestamp():
    import time
    return time.strftime("%Y%m%d-%H%M%S")
