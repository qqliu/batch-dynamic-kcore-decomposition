
with open('latencies_dblp') as f:
    lines = f.readlines()
    x = [i for i in range(counter)]
    y = [float(line) for line in lines]

fig = plt.figure()

ax1 = fig.add_subplot(111)
ax1.set_title("DBLP latencies plot")
ax1.set_xlabel("point")
ax1.set_ylabel("latency")

ax1.plot(x, y, c='r', label='data')

leg = ax1.legend()
