import re
import subprocess

int_re = r"\d+"
float_re = r"(?:\d+\.\d+|\d+\.\d+e[+-]?\d+)"
percent_re = r"\d+\.\d+%"
space_re = r"\s+"

class Event:
    def __init__(self, name, count, total_time):
        self.name = name
        self.count = count
        self.total_time = total_time

    def __repr__(self):
        return f"Event(name={self.name}, count={self.count}, total_time={self.total_time})"

    @classmethod
    def from_string(cls, event_str):
        match = re.match(
            f"(.*?){space_re}({int_re}){space_re}{float_re}{space_re}({float_re}){space_re}{float_re}{space_re}{float_re}{space_re}{float_re}{space_re}{float_re}{space_re}{float_re}{space_re}{float_re}.*",
            event_str
        )

        name = match.group(1)
        count = int(match.group(2))
        total_time = float(match.group(3))
        return cls(name, count, total_time)


# NOTE: This would ideally parse the entire file and give it back.
def profile_script(cmd, stage="Main Stage"):
    if isinstance(cmd, str):
        cmd = cmd.split()
    else:
        cmd = list(cmd)

    assert "-log_view" not in cmd

    cmd.append("-log_view",)
    result = subprocess.run(cmd, capture_output=True, text=True)

    events = collect_stage_events(result.stdout, stage)
    return events


def collect_stage_events(log_data: str, stage_name: str):
    events = []

    log_iter = iter(log_data.splitlines())
    while True:
        log_line = next(log_iter)
        if re.match(f"--- Event Stage \\d+: {stage_name}", log_line) is not None:
            for log_line in log_iter:
                if log_line.startswith("---"):
                    return tuple(events)

                elif log_line and not log_line.isspace():
                    event = Event.from_string(log_line)
                    events.append(event)
