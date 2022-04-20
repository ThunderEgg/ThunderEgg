#!env python3
import sys
import os

for file in sys.stdin:
    os.system("/Users/scott/Developer/sftwr/llvm-git/bin/clang-format -i " + file)
    f = open(file[0:-1], "r+")
    lines = f.readlines()
    f.seek(0)
    insert_i = 0
    for i in range(len(lines)):
        line = lines[i]
        if line == "{\n":
            insert_i = i+1
        if "GENERATE" in line:
            line = "for(" + line.replace("= GENERATE(",
                                         ":{").replace(");", "}){").replace("as<std::string>{},", "")
            lines.pop(i)
            lines.insert(insert_i, line)
            insert_i += 1
        if line == "}\n":
            lines[i] = "}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}\n"

    if lines[-1] == "}":
        lines[-1] = "}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}"

    for line in lines:
        f.write(line)
    f.truncate()
    f.close()
    os.system("/Users/scott/Developer/sftwr/llvm-git/bin/clang-format -i " + file)
    f = open(file[0:-1], "r+")
    lines = f.readlines()
    f.seek(0)
    prev_line = False
    lines[-1] = lines[-1] + "\n"
    for i in range(len(lines)):
        curr_line = lines[i] == "}\n"
        if curr_line and prev_line:
            lines[i] = ""
        prev_line = curr_line

    for line in lines:
        f.write(line)
    f.truncate()
    f.close()
    os.system("/Users/scott/Developer/sftwr/llvm-git/bin/clang-format -i " + file)
