import pandas as pd


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


def parse_csv(csv_file, bandwidth):
    batch_info = pd.read_csv(csv_file)
    total_time = batch_info["total_time"]
    ncalls = batch_info["count"]
    average_time = total_time / ncalls
    nflops = batch_info["nflops"]
    pessimal_memory = batch_info["pessimal_memory"]

    batch_info["flop_rate"] = flop_rate(nflops, average_time)

    batch_info["optimal_ai"] = arithmetic_intensity(nflops, batch_info["optimal_memory"])
    batch_info["optimal_peak_flops"] = peak_flop_rate(batch_info["optimal_ai"], bandwidth)

    batch_info["pessimal_ai"] = arithmetic_intensity(nflops, pessimal_memory)
    batch_info["pessimal_peak_flops"] = peak_flop_rate(batch_info["pessimal_ai"], bandwidth)
    return batch_info
