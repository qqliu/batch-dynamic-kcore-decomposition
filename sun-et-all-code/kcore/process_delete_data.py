batch_size = 100
filename = "livejournal"
out = "hua-workers-out"

files = ["../../" + out + "/" + filename + "+0_60_" +
        str(batch_size), "../../" + out + "/" + filename + "+1_60_" + str(batch_size)]
        #str(batch_size),"../../hua-out/" + filename + "+2_60_" + str(batch_size)]

best_total = 0
best_max = 0
best_average = -1
best_correct = 0
best_max_correct = 0
best_avg_correct = 0
best_batch = 0
best_avg_batch = 0
best_total_batch = 0

for f in files:
    total = 0
    total_correct = 0
    num_iterations = 0
    max_correct = 0
    max_time = 0
    max_batch = 0
    batch_avg = 0
    batch_total = 0
    with open(f, "r") as data:
        for line in data:
            time = line.split()
            if "Correct" in time:
                total_correct += float(time[2])
                if float(time[2]) > max_correct:
                    max_correct = float(time[2])
            elif "Batch" in time:
                batch_total += float(time[2])
                if float(time[2]) > max_batch:
                    max_batch = float(time[2])
            else:
                total += float(time[1])
                if float(time[1]) > max_time:
                    max_time = float(time[1])
            num_iterations += float(1.0/3)
    average_runtime = float(total)/num_iterations
    average_correct = float(total_correct)/num_iterations
    batch_avg = float(batch_total)/num_iterations
    if average_runtime < best_average or best_average == -1:
        best_average = average_runtime
        best_max = max_time
        best_total = total
        best_max_correct = max_correct
        best_correct = total_correct
        best_avg_correct = average_correct
        best_batch = max_batch
        best_avg_batch = batch_avg
        best_total_batch = batch_total

print "Batch size: " + str(batch_size)
print "Best total: " + str(best_total/1000)
print "Best avg: " + str(best_average/1000)
print "Best max: " + str(best_max/1000)
print "Best correct total: " + str(best_correct/1000)
print "Best avg correct: " + str(best_avg_correct/1000)
print "Best max correct: " + str(best_max_correct/1000)
print "Best batch total: " + str(best_total_batch/1000)
print "Best avg batch: " + str(best_avg_batch/1000)
print "Best max batch: " + str(best_batch/1000)
