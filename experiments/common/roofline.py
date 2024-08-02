def arithmetic_intensity(nflops, memory):
    return nflops / memory


def flop_rate(nflops, time):
    return nflops / time


def peak_flop_rate(arithmetic_intensity, bandwidth):
    """
    arithmetic_intensity: FLOPs/B
    bandwidth: MB/s
    """
    return arithmetic_intensity * bandwidth * 1e6
