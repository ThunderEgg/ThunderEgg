#!env python3
import sys
import os

for file in sys.stdin:
    os.system("/Users/scott/Developer/sftwr/llvm-git/bin/clang-format -i " + file)
    f = open(file[0:-1], "r+")
    lines = f.readlines()
    f.seek(0)
    for i in range(len(lines)):
        line = lines[i]
        level = ""
        if "REQUIRE" in line:
            level = "REQUIRE"
        if "CHECK" in line:
            level = "CHECK"
        if "WARN" in line:
            level = "WARN"
        if level != "":
            if " "+level+"(" in line:
                if " == " in line:
                    lines[i] = line.replace(
                        level, level+"_EQ").replace(" == ", ",")
                elif " != " in line:
                    lines[i] = line.replace(
                        level, level+"_NE").replace(" != ", ",")
                elif " > " in line:
                    lines[i] = line.replace(
                        level, level+"_GT").replace(" > ", ",")
                elif " < " in line:
                    lines[i] = line.replace(
                        level, level+"_LT").replace(" < ", ",")
                elif " >= " in line:
                    lines[i] = line.replace(
                        level, level+"_GE").replace(" >= ", ",")
                elif " <= " in line:
                    lines[i] = line.replace(
                        level, level+"_LE").replace(" <= ", ",")
                else:
                    lines[i] = line.replace(level, level+"_UNARY")
            elif " "+level+"_FALSE(" in line:
                lines[i] = line.replace(level, level+"_UNARY")

    for line in lines:
        f.write(line)
    f.truncate()
    f.close()
    os.system("/Users/scott/Developer/sftwr/llvm-git/bin/clang-format -i " + file)
