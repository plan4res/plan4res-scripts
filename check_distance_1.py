import sys
import pandas as pd
import numpy as np

# computes the distance between 1 and the solution of investment_solver
# csv_file : results_invest/Solution_OUT.csv
def check_distances(csv_file, epsilon):
    # Load the CSV file into a DataFrame
    data = pd.read_csv(csv_file, header=None)
    # Retrieve the column of numbers (the first column) in the CSV
    numbers = data[data.columns[0]]
    # Create a vector of ones with the same length as the numbers column
    ones_vector = np.ones(len(numbers))
    # Calculate the Euclidean distance between the numbers column and the ones vector
    distance = np.linalg.norm(numbers - ones_vector)

    # Check if the distance is less than the provided epsilon
    if distance < epsilon:
        print("Distance {} is less than epsilon {}".format(distance, epsilon))
        sys.exit(0)  # Exit with status 0 (success)
    else:
        print("Distance {} is greater than or equal to epsilon {}".format(distance, epsilon))
        sys.exit(1)  # Exit with status 1 (indicating distance >= epsilon)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python check_distances.py <path_to_csv> <epsilon>")
        sys.exit(1)  # Exit with error if the number of arguments is incorrect

    csv_file = sys.argv[1]
    epsilon = float(sys.argv[2])

    check_distances(csv_file, epsilon)
