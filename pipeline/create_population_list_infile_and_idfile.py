#!/usr/bin/python

import sys

# get donors and recipients from commandline
def parse_commandline (cmd_line):
  donors = []
  recipients = []
  check = True
  
  i = 0
  while i < len(cmd_line):
    if cmd_line[i] == "-d":
      while i + 1 < len(cmd_line) and cmd_line[i + 1][0] != "-":
        donors.append (cmd_line[i + 1])
        i += 1
    elif cmd_line[i] == "-r":
      while i + 1 < len(cmd_line) and cmd_line[i + 1][0] != "-":
        recipients.append (cmd_line[i + 1])
        i += 1
    else:
        sys.stderr.write ("ERROR: Unknown commandline parameter: " + cmd_line[i] + "\n")                        
    i += 1  
  
  if len(donors) == 0: 
    sys.stderr.write ("ERROR: No donors specified.\n")
    check = False
  if len(recipients) == 0:
    sys.stderr.write ("ERROR: No recipients specified.\n")
    check = False
      
  return (check, donors, recipients)

# creates populations list file for chromopainter
def create_population_list (prefix, donors, recipients):
  f = open (prefix + ".poplist", "w")
  for d in donors: f.write (d + "\tD\n")
  for r in recipients: f.write (r + "\tR\n")
  f.close()

# count populations
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

# creates idfile for chromopainter
def create_id_file (prefix, donors, recipients):
  fw = open (prefix + ".idfile", "w") 
  fr = open (prefix + ".fam", "r")
  stat = {}
  n = 0
  for line in fr:
    ind = line.strip().split()
    if ind[0] in stat: 
      stat[ind[0]] += 1
    else:
      stat[ind[0]] = 1
    if ind[0] in donors or ind[0] in recipients:
      fw.write (ind[0] + str(stat[ind[0]]) + "\t" + ind[0] + "\t1\n")
      n += 1
    else:
      fw.write (ind[0] + str(stat[ind[0]]) + "\t" + ind[0] + "\t0\n")
  fw.close()
  fr.close()
  return n

# creates parameter file for globetrotter
def create_parameter_file (prefix, donors, recipients):
  f = open (prefix + ".param", "w")
  f.write ("prop.ind: 1\n")
  f.write ("bootstrap.date.ind: 0\n")
  f.write ("null.ind: 0\n")
  f.write ("input.file.ids: " + prefix + ".idfile\n")
  f.write ("input.file.copyvectors: " + prefix + ".chunklengths.out\n")
  f.write ("save.file.main: " + prefix + ".globetrotter.main\n")
  f.write ("save.file.bootstraps: " + prefix + ".globetrotter.boot\n")
  f.write ("copyvector.popnames:")
  for d in donors:
    f.write (" " + d)
  f.write ("\n")
  f.write ("surrogate.popnames:")
  for d in donors:
    f.write (" " + d)
  f.write ("\n")
  f.write ("target.popname: " + recipients[0] + "\n")
  f.write ("num.mixing.iterations: 5\n")
  f.write ("props.cutoff: 0.001\n")
  f.write ("bootstrap.num: 20\n")
  f.write ("num.admixdates.bootstrap: 1\n")
  f.write ("num.surrogatepops.perplot: 3\n")
  f.write ("curve.range: 1 50\n")
  f.write ("bin.width: 0.1\n")
  f.write ("xlim.plot: 0 50\n")
  f.write ("prop.continue.ind: 0\n")
  f.write ("haploid.ind: 0\n")
  f.close()

prefix = sys.argv[1]
cmd_line = sys.argv[2:]

# work-flow
r = parse_commandline (cmd_line)
n = 0
if r[0]: 
  create_population_list (prefix, r[1], r[2])
  n = create_id_file (prefix, r[1], r[2])
  create_parameter_file (prefix, r[1], r[2])
print (n)
  
