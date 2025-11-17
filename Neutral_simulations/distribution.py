import sys
import numpy as np

def run(infiles):
    vals = list()
    for file in infiles:
        # Load the numpy array from the .npy file
        data = np.load(file)
        
        # Ensure that the data is a 1D array (you can adjust this if necessary)
        if data.ndim != 1:
            print(f"Warning: {file} is not a 1D array. Skipping...", file=sys.stderr)
            continue
        
        # Drop NaN values (in place)
        data = data[~np.isnan(data)]
        data = data[np.isfinite(data)]
        
        # Append to list for concatenation later
        vals.append(data)
    
    # Concatenate all the data into one long array
    dist = np.concatenate(vals)
    
    # Calculate histogram (using 1 million bins as an example)
    nbins = int(1e6)
    hist, bin_edges = np.histogram(dist, bins=nbins)
    
    # Output histogram data to stdout
    for init_bin, end_bin, count in zip(bin_edges[:-1], bin_edges[1:], hist):
        print(f"{init_bin},{end_bin},{count}")
    
    # Print quantiles to stderr
    for q in [0.5, 0.95, 0.99, 0.999, 0.9999]:
        quantile = np.quantile(dist, q)
        print(f"Quantile {q}: {quantile}", file=sys.stderr)
    
    # Calculate and print the maximum value
    max_value = np.max(dist)
    print(f"Max: {max_value}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print(
            "Usage: python distribution.py ... > <output>",
            file=sys.stderr,
        )
        sys.exit(1)
    infiles = sys.argv[1:]
    run(infiles)
