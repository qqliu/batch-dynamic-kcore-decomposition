import os
import sys
import signal
import time
import subprocess

def signal_handler(signal,frame):
  print("bye")
  sys.exit(0)
signal.signal(signal.SIGINT,signal_handler)

def shellGetOutput(str) :
  process = subprocess.Popen(str,shell=True,stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
  output, err = process.communicate()

  if (len(err) > 0):
    raise NameError(str+"\n"+output.decode('utf-8')+err.decode('utf-8'))
  return output.decode('utf-8')

def appendToFile(out, filename):
  with open(filename, "a+") as out_file:
    out_file.writelines(out)

def benchmarkToProgramPath(benchmark):
  benchmark_dict = {
    "CLDS" : "EdgeOrientation/ConcurrentPLDS/ConcurrentLDS/LDS",
    "Delayed" : "EdgeOrientation/ConcurrentPLDS/SynchronizedReads/LDS",
    "Nonlin" : "EdgeOrientation/ConcurrentPLDS/NonlinearizableReads/LDS"
  }
  return benchmark_dict.get(benchmark)

def benchmarkToIsDynamic(benchmark):
  benchmark_dict = {
    "CLDS" : True,
    "Delayed" : True,
    "Nonlin" : True
  }
  return benchmark_dict.get(benchmark)

def main():
  git_init_process = subprocess.Popen("git init ..",shell=True,stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
  output, err = git_init_process.communicate()
  git_init_delete_process = subprocess.Popen("rm -rf ../.git/",shell=True,stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
  output, err = git_init_delete_process.communicate()
  git_init_delete_process = subprocess.Popen("git submodule update --init --recursive",shell=True,stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
  output, err = git_init_delete_process.communicate()

  # Read parameters from setup file
  with open('read.txt') as parameters_file:
    for line in parameters_file:
      line = line.strip()
      split = [x.strip() for x in line.split(':')]
      if len(split) <= 1:
        continue
      params = [x.strip() for x in split[1].split(',')]
      if line.startswith("Input graph directory"):
        read_dir = split[1]
      elif line.startswith("Dynamic graphs"):
        files = params.copy()
      elif line.startswith("Output directory"):
        write_dir = split[1]
      elif line.startswith("Benchmarks"):
        programs = [benchmarkToProgramPath(x) for x in params]
        is_dynamic = [benchmarkToIsDynamic(x) for x in params]
        program_pres = params.copy()
      elif line.startswith("Numbers of workers"):
        num_workers = params.copy()
      elif line.startswith("Deltas"):
        epss = params.copy()
      elif line.startswith("Percentile"):
        percentiles = params.copy()
      elif line.startswith("Initial Graph File"):
        initial_graphs = params.copy()
      elif line.startswith("Lambdas"):
        deltas = params.copy()
      elif line.startswith("Number of Levels Divisor"):
        divisors = params.copy()
      elif line.startswith("Batch sizes"):
        batch_sizes = params.copy()
      elif line.startswith("Output stats"):
        if split[1] == "True":
          stats = "-stats"
        else:
          stats = ""
      elif line.startswith("Nonlinearizable"):
        if split[1] == "True":
          nonlin = "-nonlin"
        else:
          nonlin = ""
      elif line.startswith("Num Reader Threads"):
          num_reader_threads = params.copy()
      elif line.startswith("Output sizes"):
        if split[1] == "True":
          size = "-size"
        else:
          size = ""
      elif line.startswith("Opt"):
        if split[1] == "True":
          opt = "-ins-opt"
        else:
          opt = ""
  # Setup other parameters
  program_dir = "../../benchmarks/"
  empty = "../../benchmarks/EdgeOrientation/ConcurrentPLDS/ConcurrentLDS/empty_h"
  for program in programs:
    program_path = os.path.join(program_dir, program)
    program_local_dir = os.path.dirname(program_path)
  for file_idx, filename in enumerate(files):
    for program_idx, program in enumerate(programs):
      for e in epss:
        for d in deltas:
          for divisor in divisors:
            for b in batch_sizes:
              for nw in num_workers:
                  for nt in num_reader_threads:
                      for percent in percentiles:
                        num_rounds = 3
                        out_path_components = [program_pres[program_idx], filename, e,
                            d, b, nw, divisor, nt, percent, nonlin, stats, ".out"]
                        out_filename = os.path.join(write_dir, "_".join(out_path_components))
                        program_path = os.path.join(program_dir, program)
                        ss = ("PARLAY_NUM_THREADS=" + str(nw) + " " + program_path + " "
                        "-s -i " + read_dir + filename + " -eps " + e + " "
                        "-delta " + d + " -b " + b + " " + "-readers " + nt + " "
                        + stats + " " + nonlin + " " + size + " "
                        + "-percentile " + str(percent) + " " +
                        opt + " -opt "  + str(divisor) + " " +
                        "-rounds " + str(num_rounds))
                        if len(initial_graphs) > file_idx and len(initial_graphs[file_idx]) > 0:
                            ss += " -init_graph_file " + read_dir + initial_graphs[file_idx]
                        ss += " " + empty
                        out = shellGetOutput(ss)
                        appendToFile(out, out_filename)

if __name__ == "__main__":
  main()
