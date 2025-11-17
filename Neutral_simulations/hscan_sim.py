#!/bin/env python
import msprime
import argparse, tempfile, os, subprocess, yaml
import pandas as pd
import numpy as np

def generate_demography(params):
    """
    Define demography for following simulations
    """
    demography = msprime.Demography()
    demography.add_population(name="Ancient", initial_size=params["Ne_Anc"])
    demography.add_population(name="A", initial_size=params["Ne_A"])
    demography.add_population(name="B", initial_size=params["Ne_B"])
    demography.add_population(name="AB", initial_size=params["Ne_AB"])
    demography.add_population(name="C", initial_size=params["Ne_C"])
    demography.add_population_split(time=params["T_AB"], derived=["A", "B"], ancestral="AB")
    demography.add_population_split(time=params["T_ABC"], derived=["AB", "C"], ancestral="Ancient")
    return demography

def sim_windows(demography, params, seed):
    """
    Simulate one window. This function returns an iterator. 
    """
    samples=[
        # Sample from A, B and C in the present
        msprime.SampleSet(params["Num_A"], population="A", time=0),
        msprime.SampleSet(params["Num_B"], population="B", time=0),
        msprime.SampleSet(params["Num_C"], population="C", time=0)
    ]
    
    # Compute how many windows we need to simulate an entire genome
    n_windows = int(params["L"] // params["window_size"])
    
    
    replicates = msprime.sim_ancestry(
        demography=demography,
        samples = samples,
        sequence_length=params["window_size"],
        recombination_rate=params["RHO"],
        random_seed=seed, 
        num_replicates=n_windows
    )
    for i, ts in enumerate(replicates):
        print(f"Iteration: {i}")
        yield msprime.sim_mutations(ts, rate=params["MU"], random_seed=seed)

def compute_hscan(ts):
    """
    Compute the HSCAN values as done with real data but from a tskit object
    """
    # Define some helper functions
    def ind_name(ts, ind):
        pop = ts.population(ind.population).metadata["name"]
        return f"ts_{pop}_{ind.id}"
    # Define new names 
    individual_names = [ind_name(ts, ind) for ind in ts.individuals()]
    # Write temporary VCF and sample names
    with tempfile.TemporaryDirectory() as tempdir:
        vcf_path = os.path.join(tempdir, "out.vcf")
        with open(vcf_path, 'w') as vcf:
            ts.write_vcf(vcf, individual_names=individual_names, allow_position_zero=True)
        sample_names_path = os.path.join(tempdir, "sample_names.txt")
        with open(sample_names_path, "w") as file:
            file.writelines("\n".join(individual_names))
        # Write samples from each population into a different file
        # Assuming populations are correctly named
        pop_files = {}
        populations = ["A", "B", "C"]
        for pop in populations:
            pop_file = os.path.join(tempdir, f"{pop}.txt")
            pop_files[pop] = pop_file
            subprocess.run(
                ["grep", f'ts_{pop}_', sample_names_path],
                stdout=open(pop_file, "w"),
                check=True
            )
        # Compute Fst values with vcftools
        hscan_input = os.path.join(tempdir, "output.tsv")
        bcftools_cmd = f'''bcftools view -Ov {vcf_path} -S {pop_files["A"]} | bcftools query -f '%POS[,%TGT]\n' | awk '{{gsub(/\\./,"N")}}1' |  awk '{{gsub(/\\|/,",")}}1' > {hscan_input}'''
        os.system(bcftools_cmd)
        if not os.path.exists(hscan_input):
            raise FileNotFoundError(f"bcftools command did not produce the expected output file: {hscan_input}, {bcftools_cmd}")
        out_hscan = os.path.join(tempdir, "H_out.tsv")
        hscan_cmd = f"/data/antwerpen/grp/asvardal/software/H-scan/H-scan -i {hscan_input} -d 1 > {out_hscan}"
        os.system(hscan_cmd)
        if not os.path.exists(out_hscan):
            raise FileNotFoundError(f"H-scan command did not produce the expected output file: {out_hscan}, {hscan_cmd}")
        df2 = pd.read_csv(out_hscan, delimiter="\t") 
        return pd.DataFrame(df2['H'])

    
def main(config_file, seed, outfile):
    # Read configuration file and parse seed
    with open(config_file) as stream:
        try:
            parameters = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    seed = int(seed)
    # Define demography object
    demography = generate_demography(parameters)
    # Simulate independent windows
    replicates = sim_windows(demography, parameters, seed)
    # Simulations run sequentially:
    df = pd.concat([compute_hscan(mts) for mts in replicates])
    # Drop undefined PBS values
    df = df[df["H"]>=0.0]
    df = df.dropna()
    # Save into binary format
    with open(outfile, 'wb') as f:
        np.save(f, df["H"].to_numpy())
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate PBS values under a null-model")
    parser.add_argument("config_file", help="Path to the YAML configuration file.")
    parser.add_argument("seed", type=int, help="Random seed for the simulation.")
    parser.add_argument("outfile", help="Output file for the HSCAN values.")
    
    args = parser.parse_args()
    
    # Run the main function with parsed arguments
    main(args.config_file, args.seed, args.outfile)