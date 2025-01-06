import pandas as pd
import numpy as np
import os

# Function to check distances between initial capacities and invested capacities
# csv_directory: repo csv_invest
# solution_out_file: results_invest/Solution_OUT.csv
# epsilon: convergence limit
def check_distances(csv_directory, solution_out_file, epsilon):
    try:
        # Load the solution output file into a DataFrame
        sol_invest = pd.read_csv(solution_out_file, header=None)
    except FileNotFoundError:
        print(f'Error: File {solution_out_file} not found.')
        return

    # Define file paths
    TU_file = os.path.join(csv_directory, 'TU_ThermalUnits.csv')
    RES_file = os.path.join(csv_directory, 'RES_RenewableUnits.csv')
    STS_file = os.path.join(csv_directory, 'STS_ShortTermStorage.csv')
    IN_file = os.path.join(csv_directory, 'IN_Interconnections.csv')

    # Try loading all the required CSV files
    try:
        TU = pd.read_csv(TU_file)
        RES = pd.read_csv(RES_file)
        STS = pd.read_csv(STS_file)
        IN = pd.read_csv(IN_file)
    except FileNotFoundError as e:
        print(f'Error loading file: {e}')
        return

    capacity_before = []
    capacity_after = []
    index_sol_invest = 0  # Local variable to track index in sol_invest

    # Function to calculate capacities
    def calculate_capacities(df, column_name, sol_invest, index_sol_invest):
        for _, row in df.iterrows():
            if row['MaxAddedCapacity'] > 0 or row['MaxRetCapacity'] > 0:
                capacity_before.append(row[column_name])
                try:
                    # Check if the index exists in sol_invest before accessing it
                    capacity_after_value = sol_invest.iloc[index_sol_invest, 0] * row[column_name]
                    capacity_after.append(capacity_after_value)
                    index_sol_invest += 1
                except IndexError:
                    print(f"Error: index {index_sol_invest} out of bounds for sol_invest.")
                    return index_sol_invest  # Early return if there's an error
        return index_sol_invest

    # Calculate capacities for each unit type
    index_sol_invest = calculate_capacities(TU, 'MaxPower', sol_invest, index_sol_invest)
    index_sol_invest = calculate_capacities(RES, 'MaxPower', sol_invest, index_sol_invest)
    index_sol_invest = calculate_capacities(STS, 'MaxPower', sol_invest, index_sol_invest)
    index_sol_invest = calculate_capacities(IN, 'MaxPowerFlow', sol_invest, index_sol_invest)

    # Ensure there's no division by zero
    if np.sum(capacity_before) == 0:
        print('Error: Sum of capacity_before is zero, cannot compute relative distance.')
        return

    # Calculate the relative distance
    distance = np.sum(np.abs(np.array(capacity_before) - np.array(capacity_after))) / np.sum(capacity_before)

    # Check if the distance is less than the provided epsilon
    if distance < epsilon:
        print(f'Distance {distance} is less than epsilon {epsilon}')
    else:
        print(f'Distance {distance} is greater than or equal to epsilon {epsilon}')


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print("Usage: python check_distances.py <path_to_csv_directory> <solution_out_file> <epsilon>")
        sys.exit(1)

    csv_directory = sys.argv[1]
    solution_out_file = sys.argv[2]
    epsilon = float(sys.argv[3])

    check_distances(csv_directory, solution_out_file, epsilon)
