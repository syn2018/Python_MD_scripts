import os
import urllib

if os.name == 'posix':
  os.chdir("/tmp")
  try:
     os.mkdir("PythonSoftware")
  except:
     pass	# Directory already there?
  os.chdir("PythonSoftware")
  urllib.urlretrieve("http://www.physics.cornell.edu/sethna/StatMech/ComputerExercises/PythonSoftware.tgz", "PythonSoftware.tgz")
  print "Untarring"
  os.system("tar -xvzf PythonSoftware.tgz")
  print "Moving"
  os.system("mv */* .")
  os.system("xterm -e python DoPythonDemos.py")

elif os.name == 'nt':
  try:
     os.mkdir("PythonSoftware")
  except:
     pass	# Directory already there?
  os.chdir("PythonSoftware")
  urllib.urlretrieve("http://www.physics.cornell.edu/sethna/StatMech/ComputerExercises/PythonSoftware.zip", "PythonSoftware.zip")
  print "Unzipping"
  os.system("Unzip it")
  print "Moving"
  os.system("move */* .")
