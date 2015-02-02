# This scrips extracts timings from the output of a gdem run

import re
import sys

total_times = dict()
comm_times = dict()

total_time_re = re.compile("\\(Process (\d)\\) Elapsed time: (\d+\\.\d+)")
comm_time_re = re.compile("\\(Process (\d)\\) Time .*: (\d+\\.\d+)")
diam_re = re.compile("Effective diameter .*")

def process_line(line):
    m = total_time_re.match(line)
    if m:
        total_times[int(m.group(1))] = float(m.group(2))
        return
    m = comm_time_re.match(line)
    if m:
        key = int(m.group(1))
        if key not in comm_times:
            comm_times[int(m.group(1))] = float(m.group(2))
        else:
            comm_times[int(m.group(1))] += float(m.group(2))


if __name__ == "__main__":
    for line in sys.stdin:
        if diam_re.match(line):
            print line
        process_line(line)

    elapsed = max(total_times.values())
    comm = max(comm_times.values())
    perc = comm / elapsed * 100

    print "Total elapsed time:", elapsed, "seconds"
    print "Communication time:", comm, "seconds"
    print "Percentage of communication:", perc, "%"

