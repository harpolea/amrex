"""
Converts the addresses in the backtrace to file names and line numbers.
"""

import re
import subprocess
import sys

if __name__ == "__main__":
    if len(sys.argv) > 1:
        input_filename = sys.argv[1]
    else:
        input_filename = "Backtrace.rg_0_rl_0"

    f = open(input_filename, "r+")

    f_content = f.readlines()

    executable = None


    for i, line in enumerate(f_content):
        if executable is None:
            m = re.search("\./(\S+.ex)", line)
            if m is not None:
                executable = m.group(1)
        m = re.search("\((\+0x\w+)\)", line)
        if m is not None:
            addr = m.group(1)
            bash_command = "addr2line -Cfie " + executable + ' ' + addr
            process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()
            process.terminate()
            output = output.decode('utf-8').split('\n') # convert bytes to string
            f_content[i+1] = '\t' + output[0]
            f_content[i+2] = '\n\t' + output[1]

    print(f_content)
    f.seek(0)
    f.truncate()
    f.write('\n'.join(f_content))
