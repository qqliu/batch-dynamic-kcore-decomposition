import os
import sys
import signal
import time
import subprocess

def signal_handler(signal,frame):
  print "bye\n"
  sys.exit(0)
signal.signal(signal.SIGINT,signal_handler)

def shellGetOutput(str) :
  process = subprocess.Popen(str,shell=True,stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
  output, err = process.communicate()


  #if (len(err) > 0):
  #    raise NameError(str+"\n"+output+err)
  return output

def computeTimeout(out):
  time = 0
  for line in out.splitlines():
    line = line.strip()
    split = [x.strip() for x in line.split(':')]
    if split[0].startswith("### Batch Running Time"):
      time += float(split[1])
  return time

def appendToFile(out, filename):
  with open(filename, "a+") as out_file:
    out_file.writelines(out)

# cluster time, total time, modularity, precision, recall, # comm,
# min comm, max comm, avg comm
def main():
  # Configured for Test 1
<<<<<<< HEAD
  files = ["livejournal"] #,"livejournal_insertion_edges","youtube_insertion_edges","orkut_insertion_edges"]
  batch_sizes = [10000]#, 1000, 10000, 100000, 1000000, 10000000]
  num_workers = [60]#[1, 2, 4, 8, 16, 30, 60]
=======
  files = ["youtube", "orkut"] #,"livejournal_insertion_edges","youtube_insertion_edges","orkut_insertion_edges"]
  batch_sizes = [1000000]
  num_workers = [60]
  #num_workers = [1, 2, 4, 8, 16, 30, 60]
>>>>>>> 40d261e2d4c91049da3c85fa2fa4b187e4b8b6ad
  read_dir = "/home/qliu19/dynamic_graph/"
  write_dir = "/home/qliu19/hua-big-out/"
  num_rounds = 3
  for threads in num_workers:
    for f in files:
        for batch_size in batch_sizes:
            for r in range(num_rounds):
                ss = "./kcore " + read_dir + "hua_" + f + "_insertion_edges" + " " + str(threads) + " " + read_dir + "hua_" + f + "_deletion_edges " + str(batch_size)
                out_filename = write_dir + f + "+" + str(r) + "_" + str(threads) + "_" + str(batch_size)
                print ss
                out = shellGetOutput(ss)
                #print out
                appendToFile(out, out_filename)

if __name__ == "__main__":
  main()
