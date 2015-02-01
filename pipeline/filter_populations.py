#!/usr/bin/python

import sys

## functions
def collect_statistics (file, populations):
  stats = {}
  used = []
  if len(populations) == 0: use_all_populations = True
  else: use_all_populations = False
  
  for line in file:
    line = line.strip()
    pop = line.split('\t')[0]
    if not use_all_populations and pop not in populations: continue
    if pop in stats: 
      stats[pop] += 1
    else:
      used.append(pop) 
      stats[pop] = 1 
      
  for p in populations:
    if p not in used: sys.stderr.write ("WARNING: Population " + p + " not in the data set!\n")  
  return stats
  
def print_stats (stats):
    for key, value in stats.items():
      print(key + "\t" + str(value))
    return


## get commandline parameters
INPUT = sys.argv[1]
populations = []
get_statistics = False

for i in range (2, len(sys.argv)):
  if sys.argv[i] != "--stat":
    populations.append(sys.argv[i])
  else:
    get_statistics = True


f = open(INPUT, "r")

## main workflow
# getting population statistics only
if get_statistics:
  stats = collect_statistics (f, populations)
  print_stats (stats)
# filtering populations workflow  
elif len(populations) != 0:
  used = []
  for line in f:
    line = line.strip()
    pop = line.split('\t')[0]
    if pop in populations:
      if pop not in used: used.append(pop)
      print(line)
      
  for p in populations:
    if p not in used: sys.stderr.write ("WARNING: Population " + p + " not in the data set!\n")

f.close()