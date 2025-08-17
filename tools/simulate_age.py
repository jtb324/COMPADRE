import pandas as pd
import numpy as np
import sys

"""

This script takes a positional .fam file argument as input and returns 

"""

def read_fam_file(file_path):
    columns = ['FID', 'IID', 'Father', 'Mother', 'Sex', 'Phenotype']
    return pd.read_csv(file_path, sep='\s+', header=None, names=columns)

def assign_initial_generation(pedigree):
    generations = {}
    
    def assign_gen(iid, gen=0):
        if iid in generations:
            return
        generations[iid] = gen
        children = pedigree[(pedigree['Father'] == iid) | (pedigree['Mother'] == iid)]['IID']
        for child in children:
            assign_gen(child, gen + 1)
    
    founders = pedigree[(pedigree['Father'] == '0') & (pedigree['Mother'] == '0')]['IID']
    for founder in founders:
        assign_gen(founder)
    
    return generations

def refine_generations(pedigree, initial_generations):
    refined_generations = initial_generations.copy()
    
    gen_to_iids = {}
    for iid, gen in refined_generations.items():
        gen_to_iids.setdefault(gen, set()).add(iid)
    
    max_gen = max(gen_to_iids.keys())
    
    for gen in range(max_gen, -1, -1):
        for iid in gen_to_iids.get(gen, set()):
            father = pedigree.loc[pedigree['IID'] == iid, 'Father'].iloc[0]
            mother = pedigree.loc[pedigree['IID'] == iid, 'Mother'].iloc[0]
            
            if father != '0' and father in refined_generations:
                refined_generations[father] = gen - 1
            if mother != '0' and mother in refined_generations:
                refined_generations[mother] = gen - 1
    
    return refined_generations

def find_spouses(pedigree):
    spouses = {}
    for _, row in pedigree.iterrows():
        if row['Father'] != '0' and row['Mother'] != '0':
            spouses[row['Father']] = row['Mother']
            spouses[row['Mother']] = row['Father']
    return spouses

def simulate_ages(pedigree, generations, spouses, max_age=100, age_spread=5, min_parent_age=14):
    max_gen = max(generations.values())
    ages = {}

    def initial_age_assignment():
        for iid, gen in generations.items():
            age = max_age - (gen * 25) + np.random.normal(0, age_spread)
            ages[iid] = max(0, min(int(age), max_age))

    def adjust_parent_child_ages():
        for _, row in pedigree.iterrows():
            if row['IID'] in ages and row['Father'] in ages and row['Mother'] in ages:
                child_age = ages[row['IID']]
                father_age = ages[row['Father']]
                mother_age = ages[row['Mother']]
                if father_age - child_age < min_parent_age:
                    ages[row['Father']] = child_age + min_parent_age
                if mother_age - child_age < min_parent_age:
                    ages[row['Mother']] = child_age + min_parent_age

    def adjust_spouse_ages():
        for iid, spouse_iid in spouses.items():
            if iid in ages and spouse_iid in ages:
                avg_age = (ages[iid] + ages[spouse_iid]) / 2
                ages[iid] = int(avg_age + np.random.randint(-2, 3))
                ages[spouse_iid] = int(avg_age + np.random.randint(-2, 3))

    def verify_ages():
        for _, row in pedigree.iterrows():
            if row['IID'] in ages and row['Father'] in ages and row['Mother'] in ages:
                child_age = ages[row['IID']]
                father_age = ages[row['Father']]
                mother_age = ages[row['Mother']]
                if father_age - child_age < min_parent_age or mother_age - child_age < min_parent_age:
                    return False
        return True

    initial_age_assignment()
    
    max_iterations = 100
    for _ in range(max_iterations):
        adjust_parent_child_ages()
        adjust_spouse_ages()
        if verify_ages():
            break
    else:
        print("Warning: Maximum iterations reached. Some age relationships may not be ideal.")

    return ages

def main(fam_file_path):
    pedigree = read_fam_file(fam_file_path)
    initial_generations = assign_initial_generation(pedigree)
    refined_generations = refine_generations(pedigree, initial_generations)
    spouses = find_spouses(pedigree)
    ages = simulate_ages(pedigree, refined_generations, spouses)
    
    pedigree['Age'] = pedigree['IID'].map(ages)
    pedigree['Generation'] = pedigree['IID'].map(refined_generations)

    output_file = fam_file_path.replace('.fam', '_with_ages.fam')
    pedigree.to_csv(output_file, index=False)
    print(f"Results saved to: {output_file}")
    
    return pedigree



if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Error: No .fam file provided as a positional argument.\nRe-run like this: python3 simulate_age.py your_fam_file.fam")
        sys.exit(1)
    
    fam_file = sys.argv[1]
    main(fam_file)