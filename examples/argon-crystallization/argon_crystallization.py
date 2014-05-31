#!/usr/bin/python
import subprocess
import os
import os.path
import signal
from sys import argv
from argparse import ArgumentParser
import shutil


parser = ArgumentParser()
parser.add_argument("config_file")
parser.add_argument("--id", nargs='?', default="tmp")
args = parser.parse_args()

output_dir = os.path.abspath("/home/svenni/scratch/tmp/")

if args.id != "tmp":
    try:
        from sumatra.projects import load_project
        output_dir = os.path.join(os.path.abspath(load_project().data_store.root), args.id)
    except ImportError:
        pass
    
output_path = os.path.join(output_dir, "atoms*.bin")

current_path = os.path.dirname(os.path.realpath(__file__))

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
config_file = os.path.abspath(args.config_file)

build_path = os.path.abspath(os.path.join(current_path, "..", "..", "..", "build-argon-crystallization"))
project_path = os.path.abspath(os.path.join(current_path, "..", ".."))

print "Building in:\n", build_path

if os.path.exists(build_path):
    shutil.rmtree(build_path)

os.makedirs(build_path)

subprocess.call(["qmake", project_path, "CONFIG+=nogui", "CONFIG+=mpi", "CONFIG+=notests", "CONFIG+=noapp"], cwd=build_path)
subprocess.call(["make", "-j", "8"], cwd=build_path)

staterunner_path = os.path.join(build_path, "examples", "argon-crystallization")
lib_path = os.path.join("..", "..", "src", "libs")

env = dict(os.environ)
env['LD_LIBRARY_PATH'] = lib_path
#proc = subprocess.call(["./staterunner", states_file, output_file], cwd=staterunner_path, env=env)
run_argument = ["./argon-crystallization", config_file, output_path]
print " ".join(run_argument)
proc = subprocess.call(run_argument, cwd=staterunner_path, env=env)

print "Results saved to this directory:\n", output_path
