import os
from os.path import isfile, join
import subprocess

BIN = "./test.sh"
INPUTS_PATH = "./inputs/"
OUTPUTS_PATH = "./outputs/"
PASSED = 1
FAILED = 0

class bcolors:
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    FAIL = '\033[91m'
    OKGREEN = '\033[92m'



def run_test(input, expected_output):
    test = subprocess.Popen([BIN, input], stdout=subprocess.PIPE)
    output = test.stdout.read().decode()
    if expected_output == output:
        print(f"    {bcolors.BOLD}{bcolors.OKGREEN}PASSED{bcolors.ENDC}")
        return PASSED
    else:
        print(f"    {bcolors.BOLD}{bcolors.FAIL}FAILED{bcolors.ENDC}")
        print("Expected:")
        print(expected_output)
        print("Actual:")
        print(output)
        return FAILED
    

files = [
    file[:len(file)-3] for file in os.listdir(INPUTS_PATH) if isfile((join(INPUTS_PATH, file)))
]


num_of_passed = 0
for file in files:
    print(f"{bcolors.BOLD}Test: {bcolors.ENDC} {file}")
    with open(INPUTS_PATH+file+".in", "r") as input_f:
        with open(OUTPUTS_PATH+file+".out", "r") as expected_output_f:
            inpt = input_f.read()
            expected = expected_output_f.read()
            num_of_passed += run_test(inpt, expected)

print(f"\n{bcolors.BOLD}Result: {num_of_passed}/{len(files)} tests passed{bcolors.ENDC}")