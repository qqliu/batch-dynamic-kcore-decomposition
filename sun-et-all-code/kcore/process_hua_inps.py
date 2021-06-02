files = ["livejournal_deletion_edges", "dblp_deletion_edges", "orkut_deletion_edges", "youtube_deletion_edges"]

for f in files:
    hua = open("hua_" + f, "w")
    with open(f, "r") as inp:
        for line in inp:
            items = line.split()
            if items[0] == "-":
                hua.write(items[1] + " " + items[2] + "\n")
