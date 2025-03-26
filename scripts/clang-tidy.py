import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

# Purpose: Run clang-tidy on the source files

# run clang-tidy on the source files
# Usage: python clang-tidy.py [path to clang-tidy] [path to compile_commands.json] [paths to source files...]

def run_clang_tidy_on_file(clang_tidy_path, compile_commands_path, file_path):
    # Construct the clang-tidy command
    command = [
        clang_tidy_path,
        file_path,
        f"-p={compile_commands_path}",
        "-checks=-*,misc-include-cleaner",
        "-warnings-as-errors=*"
    ]
    #print(f"Running: {' '.join(command)}")
    # Run the clang-tidy command
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return file_path, result.returncode, result.stdout, result.stderr

def run_clang_tidy(clang_tidy_path, compile_commands_path, source_dirs):
    # Collect all .cpp and .h files
    files_to_check = []
    for source_dir in source_dirs:
        for root, _, files in os.walk(source_dir):
            for file in files:
                if file == "doctest.h":
                    continue  # Skip the doctest.h file
                if file.endswith(('.cpp', '.h')):
                    files_to_check.append(os.path.join(root, file))

    # Use ThreadPoolExecutor to run clang-tidy in parallel
    has_error = False
    with ThreadPoolExecutor() as executor:
        future_to_file = {
            executor.submit(run_clang_tidy_on_file, clang_tidy_path, compile_commands_path, file): file
            for file in files_to_check
        }

        for future in as_completed(future_to_file):
            file_path = future_to_file[future]
            try:
                file_path, returncode, stdout, stderr = future.result()
                if returncode != 0:
                    has_error = True
                    print(f"Error in file: {file_path}")
                    print(stderr)
                    print(stdout)
                else:
                    print(f"Passed: {file_path}")
            except Exception as e:
                has_error = True
                print(f"Exception occurred while processing {file_path}: {e}")

    # Exit with nonzero if any command failed
    if has_error:
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python clang-tidy.py [path to clang-tidy] [path to compile_commands.json] [paths to source files...]")
        sys.exit(1)

    clang_tidy_path = sys.argv[1]
    compile_commands_path = sys.argv[2]
    source_dirs = sys.argv[3:]

    run_clang_tidy(clang_tidy_path, compile_commands_path, source_dirs)
