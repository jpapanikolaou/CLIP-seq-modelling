import matplotlib.pyplot as plt
from scipy.stats import poisson
import numpy as np
import pandas as pd
def create_poisson_qq_plot(data, chromosome):
    # Filter data for the specific chromosome
    filtered_data = data[data['chromosome'] == chromosome]
    if filtered_data.empty:
        print(f"No data available for chromosome {chromosome}.")
        return

    # Extract the specific lambda estimate from the filtered data
    lambda_estimate = filtered_data['chrom_global_lambda'].iloc[0]
    print(f"Using lambda = {lambda_estimate} for chromosome {chromosome}")

    # Sort the observed data in 'counts'
    sorted_data = np.sort(filtered_data['counts'])

    # Extend the probability range to better visualize high counts
    probs = np.linspace(0.5 / len(sorted_data), 1 - 0.5 / len(sorted_data), len(sorted_data))
    poisson_quantiles = poisson.ppf(probs, lambda_estimate)

    # Create Q-Q plot
    plt.figure(figsize=(8, 8))
    plt.plot(poisson_quantiles, sorted_data, 'o', label='Observed vs. Theoretical Quantiles')
    plt.plot(sorted_data, sorted_data, 'r--', label='y = x (Perfect Fit)')
    plt.xlabel('Theoretical Poisson Quantiles')
    plt.ylabel('Observed Quantiles')
    plt.title(f'Poisson Q-Q Plot for {chromosome}')
    plt.legend()
    plt.show()


def add_global_p_value(data):
    """
    Updates the DataFrame to include a 'p_value' column, calculated using the Poisson survival function
    based on the counts for each entry and its corresponding lambda value from 'chrom_global_lambda'.

    Parameters:
        data (pd.DataFrame): DataFrame with columns 'counts' and 'chrom_global_lambda'.

    Returns:
        pd.DataFrame: The DataFrame with an additional 'p_value' column.
    """
    data['p_value'] = data.apply(
        lambda row: poisson.sf(row['counts'] - 1, row['chrom_global_lambda']),
        axis=1
    )
    return data


def add_chromosome_global_lambda(dataframe):
    """
    Calculates the global lambda (average counts) for each chromosome in the given DataFrame and
    adds this as a new column 'chrom_global_lambda' to the DataFrame.

    Parameters:
        dataframe (pd.DataFrame): DataFrame with at least two columns 'chromosome' and 'counts'.

    Returns:
        pd.DataFrame: Modified DataFrame with an additional column 'chrom_global_lambda'.
    """
    # Calculate lambda for each chromosome based on counts
    lambda_per_chromosome = dataframe.groupby('chromosome')['counts'].mean().reset_index()
    lambda_per_chromosome.columns = ['chromosome', 'chrom_global_lambda']  # Rename for clarity

    # Merge the lambda values back into the original dataframe
    modified_dataframe = pd.merge(dataframe, lambda_per_chromosome, on='chromosome', how='left')

    return modified_dataframe
