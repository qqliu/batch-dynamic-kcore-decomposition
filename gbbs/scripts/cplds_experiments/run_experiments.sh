#!/bin/bash

(cd ../../benchmarks/EdgeOrientation/ConcurrentPLDS/ConcurrentLDS/ && make)
(cd ../../benchmarks/EdgeOrientation/ConcurrentPLDS/NonlinearizableReads/ && make)
(cd ../../benchmarks/EdgeOrientation/ConcurrentPLDS/SynchronizedReads/ && make)

echo "Running experiments";
echo "-------------------------------";

echo "Latency experiments: average latency, 99% latency, 99.99% latency"
echo "-------------------------------";
python3 latency_test.py
python3 latency_read.py
echo "-------------------------------";
echo "Update time experiments: average update time and max update time"
echo "-------------------------------";

python3 update_time_read.py
echo "-------------------------------";
echo "Approximation experiments: average approximation factor and max approximation factors"
echo "-------------------------------";

python3 approx_test.py
python3 approx_read.py
echo "-------------------------------";
echo "Scalability of Read Throughputs"
echo "-------------------------------";
python3 read_test.py
python3 read_read.py
echo "-------------------------------";
echo "Scalability of Write Throughputs"
echo "-------------------------------";
python3 write_test.py
python3 write_read.py
echo "-------------------------------";
