import re
import os


# compute the distance between the cost of the solutions found at the last 2 iterations of investment_solver
# log_file: log file of investment_solver
# epsilon: convergence test
def check_distances(log_file, epsilon):
    try:
        with open(log_file, 'r') as file:
            content = file.read()
            
            # Find all instances of "Solution value" followed by a number in scientific or decimal notation
            matches = re.findall(r'Solution value:\s*([-+]?\d*\.\d+(e[+-]?\d+)?|\d+e[+-]?\d+)', content)
            
            if len(matches) < 2:
                print('Not enough solution values found in the log file.')
                return
            
            # Extract numbers from the regex matches
            last_value = float(matches[-1][0])
            second_last_value = float(matches[-2][0])
            
            # Calculate the absolute difference
            difference = abs(last_value - second_last_value)/abs(second_last_value)

            # Check if the difference is lower than the provided epsilon
            if difference < epsilon:
                print('The difference between the last two values ({:.5e} and {:.5e}) is less than epsilon {}: Difference = {:.5e}'.format(
                    last_value, second_last_value, epsilon, difference))
            else:
                print('The difference between the last two values ({:.5e} and {:.5e}) is greater than or equal to epsilon {}: Difference = {:.5e}'.format(
                    last_value, second_last_value, epsilon, difference))

    except FileNotFoundError:
        print('Error: File {} not found.'.format(log_file))
        return

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python check_distances.py <log_file> <epsilon>")
        sys.exit(1)

    log_file = sys.argv[1]
    epsilon = float(sys.argv[2])

    check_distances(log_file, epsilon)
