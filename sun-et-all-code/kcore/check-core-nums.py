from itertools import izip

with open("new_all_inserts") as julienne, open("jul_check") as hua:
    for x, y in izip(julienne, hua):
        x = x.split()
        y = y.split()
        if x[0] != y[0] or x[1] != y[1]:
            print "ERROR"
