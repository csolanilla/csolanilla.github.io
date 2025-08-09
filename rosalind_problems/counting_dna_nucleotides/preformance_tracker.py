import subprocess
import time

def run_script_and_time(script_path, iterations=10):
    """
    Executes a given Python script a specified number of times and
    tracks the execution time for each run.

    Args:
        script_path (str): The path to the Python script to execute.
        iterations (int): The number of times to execute the script.
    """
    execution_times = []
    print(f"Executing '{script_path}' {iterations} times...")

    for i in range(iterations):
        start_time = time.perf_counter()  # Use perf_counter for precise timing
        try:
            # Execute the script using subprocess.run
            # capture_output=True and text=True to capture stdout/stderr as text
            result = subprocess.run(['python', script_path], capture_output=True, text=True, check=True)
            end_time = time.perf_counter()
            elapsed_time = end_time - start_time
            execution_times.append(elapsed_time)
            print(f"Run {i+1}: Finished in {elapsed_time:.4f} seconds.")
            # Optional: Print script's output
            # print("Script Output:\n", result.stdout)
            # if result.stderr:
            #     print("Script Error:\n", result.stderr)
        except subprocess.CalledProcessError as e:
            print(f"Error running script (Run {i+1}): {e}")
            print("Stderr:", e.stderr)
        except FileNotFoundError:
            print(f"Error: Script '{script_path}' not found.")
            break # Exit if the script itself isn't found

    if execution_times:
        print("\n--- Summary ---")
        for i, t in enumerate(execution_times):
            print(f"Run {i+1} time: {t:.4f} seconds")
        print(f"Average execution time for **{script_path}**: {sum(execution_times) / len(execution_times):.4f} seconds")

# Call the function to run 'my_script.py' 10 times
run_script_and_time('solution.py', iterations=10)
print()
run_script_and_time('solutions_online.py', iterations=10)