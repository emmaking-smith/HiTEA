'''
Run an ANOVA-Tukey test on Homogenous Hydrogenation data.
'''
import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.formula.api import ols
import argparse
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import os
from pathlib import Path

def init_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_path', '-d', type=str, default='data/cleaned_datasets/hetero_hydrogenation.csv')
    parser.add_argument('--save_dir', '-s', type=str, default='.')
    parser.add_argument('--reaction_type', '-r', type=str, default='alkene', help = 'Which type of heterogeneous hydrogenation will be analyzed.')
    parser.add_argument('--time_data', '-t', type=int, default=None, help='Are you using time data or not. Default is no temporal analysis.')
    parser.add_argument('--remove_negative_rxns', '-n', type=bool, default=False, help='Do you wish to remove the 0% yield reactions (T/F)? Default is False, no 0% yield reactions removed.')
    return parser.parse_args()

def find_factor_count(tukey_df, factor):
    count = len(tukey_df.loc[tukey_df['group1'] == factor]) + len(tukey_df.loc[tukey_df['group2'] == factor])
    return count

def find_average_zscore_factor(df, column, factor):
    df_i = df.loc[df[column] == factor, ['z_score']]
    return float(np.mean(df_i))

def cat_lig_column(df):
    cat_lig = []
    for i in range(len(df)):
        cat = df.loc[i, 'catalyst']
        lig = df.loc[i, 'ligand']
        if pd.isna(lig) == False:
            cat_lig.append(cat + '+' + lig)
        else:
            cat_lig.append(cat + '+')
    df['cat_lig'] = cat_lig
    df = df.drop(columns=['catalyst', 'ligand'])
    return df
            
def one_hot_solvent(df):
    unique_solvents = df['Solvent_1_Name'].unique()
    one_hot_solvents = []
    for i in range(len(df)):
        one_hot_solvents.append([int(df.loc[i, 'Solvent_1_Name'] == x) for x in unique_solvents])
    one_hot_solvents = np.array(one_hot_solvents)
    for i, name in enumerate(unique_solvents):
        df['Solvent_1_Name_' + name] = one_hot_solvents[:, i].tolist()
    df = df.drop(columns=['Solvent_1_Name'])
    return df

def one_hot_reagent(df):
    unique_reagents = df['Reagent_1_Short_Hand'].unique()
    one_hot_reagents = []
    for i in range(len(df)):
        one_hot_reagents.append([int(df.loc[i, 'Reagent_1_Short_Hand'] == x) for x in unique_reagents])
    one_hot_reagents = np.array(one_hot_reagents)
    for i, name in enumerate(unique_reagents):
        df['Reagent_1_Short_Hand_' + name] = one_hot_reagents[:, i].tolist()
    df = df.drop(columns=['Reagent_1_Short_Hand'])
    return df

def one_hot_cat_lig(df):
    unique_catlig = df['cat_lig'].unique()
    one_hot_catlig = []
    for i in range(len(df)):
        one_hot_catlig.append([int(df.loc[i, 'cat_lig'] == x) for x in unique_catlig])
    one_hot_catlig = np.array(one_hot_catlig)
    for i, name in enumerate(unique_catlig):
        df['cat_lig_' + name] = one_hot_catlig[:, i].tolist()
    df = df.drop(columns=['cat_lig'])
    return df

homo_rxn_type_dict = {'alkene' : ['c12c(ncnc1N1CCN(CC1)C(OC(C)(C)C)=O)C=C(CC2)c1c(ccc2c1cnn2[C@@H]1OCCCC1)C',
                                   'C1C(C(=C1c1ccccc1)C)=O',
                                   'c1cccc(c1)C(\\C(F)(F)F)=C/C(O)=O',
                                   'c1nccnc1C(=NCc1ccc(cc1)OC)C',
                                   'n1(nc(c(c1)C=1CCCN1)C)C',
                                   'C(C(F)(F)F)(C(O)=O)=C',
                                   'c1(cnc2c(c1C)c(cn2C)C1=CCN(C(C1)(C)C)C(OC(C)(C)C)=O)[N+]([O-])=O',
                                   'c12c(C(NCC1=C)=O)cccc2',
                                   'C1C=C(CC=C1O[Si](C)(C)C(C)(C)C)C',
                                   'N1(C(=N[C@](C(S1(=O)=O)=C)(C)c1scc(n1)NC(c1ccc(cn1)OC(F)F)=O)NC(OC(C)(C)C)=O)C',
                                   'c1(cc(c(nn1)N)C(c1ccccc1)=C)Cl'],
                        
                        'CO_reduction' : ['c1(cc(c2c(c1Cl)C(N(CC2)Cc1c(cc(nc1OCc1ccccc1)C)OC)=O)C)C(=O)[C@@H]1COCC1',
                                         'c1cc(c(c(c1I)C(C)=O)Cl)F',
                                         'O1[C@H]([C@H]2[C@@H]([C@H]1C(=O)c1cc(c(cc1)F)F)OC(O2)(C)C)n1cc(c2c1ncnc2Cl)F',
                                         'c1(c(ccc(c1)F)OC)C(C)=O',
                                         'c1(cc(c2c(c1Cl)C(N(CC2)Cc1c(cc(nc1OCc1ccccc1)C)C)=O)Cl)C([C@@H]1CCOC1)=O',
                                         'c1(cc(c2c(c1Cl)C(N(CC2)Cc1c(cc(nc1OCc1ccccc1)C)C)=O)Cl)C(=O)C1CCN(CC1)C(=O)OC(C)(C)C']
                        
                       }

def filter_out_rxn_types(df, rxn_type):
    rxns = homo_rxn_type_dict[rxn_type]
    filtered_df = pd.DataFrame()
    for i in range(len(df)):
        rxn = df.loc[i, 'Reactant_1_SMILES']
        if rxn in rxns:
            filtered_df = pd.concat((filtered_df, df[i:i+1]))
    return filtered_df.reset_index(drop=True)

def zscore(df, product):
    product_rows = df.loc[df['PRODUCT_STRUCTURE'] == product]
    if len(product_rows) > 1:
        z_scores = []
        mean = np.mean(product_rows['Product_Yield_PCT_Area_UV'])
        std = np.std(product_rows['Product_Yield_PCT_Area_UV'])
        
        for i in product_rows['Product_Yield_PCT_Area_UV']:
            z_scores.append( (i - mean) / (std + 0.01) )
        df.loc[list(product_rows.index), 'z_score'] = z_scores

    else:
        df = df.drop(list(product_rows.index)).reset_index(drop=True)
    return df
    
def main():
    # Initialize the arguments.
    args = init_args()
    data_path = args.data_path
    save_dir = args.save_dir
    reaction_type = args.reaction_type
    time_data = args.time_data
    remove_negative_rxns = args.remove_negative_rxns

    Path(save_dir).mkdir(parents=True, exist_ok=True)
    
    if reaction_type not in list(homo_rxn_type_dict.keys()):
        print(f'{reaction_type} is not in listed reactions. Options are {list(homo_rxn_type_dict.keys())}. Modify the rxn_type_dict or use one of the specified reaction types.')
        assert reaction_type in list(homo_rxn_type_dict.keys())
    
    # Import data.
    df = pd.read_csv(data_path)
    df = df.loc[pd.isna(df['Reaction_T']) == False]
    df = df.loc[pd.isna(df['Reaction_Time_hrs']) == False].reset_index(drop=True)
    df.loc[pd.isna(df['Reagent_1_Short_Hand']) == True, ['Reagent_1_Short_Hand']] = 'None'
    df.loc[df['ligand'] == 'R-JosiphosSL-J009-1', ['ligand']] = '"R-Josiphos SL-J009-1'
    df = df.loc[pd.isna(df['catalyst']) == False]
    
    
    df = filter_out_rxn_types(df, reaction_type)
    if remove_negative_rxns == True:
        df = df.loc[df['Product_Yield_PCT_Area_UV'] > 0]
        
    if time_data is not None:
        df = df.loc[df['Year'] == time_data]
        
    df = df[['Product_Yield_PCT_Area_UV', 'PRODUCT_STRUCTURE', 'Solvent_1_Name', 'Reaction_T', 'Reactant_1_SMILES', 
             'catalyst', 'ligand', 'Reagent_1_Short_Hand']]
             
    df = df.reset_index(drop=True)
    
    if len(df) > 100:

        # Processing the data.
        df = cat_lig_column(df)

        # Adding the z-score info
        for product in df['PRODUCT_STRUCTURE'].unique():
            df = zscore(df, product)
        
        # ANOVA
        model = ols('z_score ~ C(cat_lig) + C(Reagent_1_Short_Hand)', data=df).fit()
        anova = sm.stats.anova_lm(model, type=2)
        statistically_significant_feats = list(anova.loc[anova['PR(>F)'] < 0.05].index)
        print(f'Statistically significant features: {statistically_significant_feats}')
        if len(statistically_significant_feats) > 0:
        # Perform Tukey tests on statistically significant feats.
            for feature in statistically_significant_feats:
                if feature == 'C(cat_lig)' in statistically_significant_feats:
                    print('Performing Tukey Test on Catalyst/Ligands...')
                    tukey_cats = pairwise_tukeyhsd(df['z_score'], df['cat_lig'])

                    # Forming the tukey dataframe.
                    tukey_df = pd.DataFrame(data=tukey_cats._results_table.data[1:], columns=tukey_cats._results_table.data[0])
                    tukey_df = tukey_df.loc[tukey_df['p-adj'] < 0.05].reset_index(drop=True)

                    if len(tukey_df) > 0:

                        # Finding top 15 outlier catalysts
                        unique_factors = np.unique(tukey_df['group1'].unique().tolist() + tukey_df['group2'].unique().tolist())
                        counts = [find_factor_count(tukey_df, x) for x in unique_factors]
                        top_15_factor_performance = pd.DataFrame({'factor' : unique_factors, 'count' : counts})
                        top_15_factor_performance = top_15_factor_performance.sort_values(by=['count'])[::-1][0:15]

                        # Finding the average z-score for each top outlier catalyst.
                        average_z_scores = [find_average_zscore_factor(df, 'cat_lig', x) for x in top_15_factor_performance['factor']]
                        top_15_factor_performance['average_z_score'] = average_z_scores
                        top_15_factor_performance = top_15_factor_performance[['factor', 'average_z_score']]

                        # Saving the dataframe.
                        if time_data is not None:
                            name = reaction_type + '_homo_hydrogenation_top_15_outlier_catalysts_' + str(time_data) + '.csv'
                        elif remove_negative_rxns == True:
                            name = reaction_type + '_homo_hydrogenation_top_15_outlier_catalysts_no_0%_yield.csv'
                        else:
                            name = reaction_type + '_homo_hydrogenation_top_15_outlier_catalysts.csv'
                        top_15_factor_performance.to_csv(os.path.join(save_dir, name))
                        print(f'Saved outlier catalyst dataframe to {os.path.join(save_dir, name)}.')

                    else:
                        print('No catalyst was statistically better than any other.')
                else:   
                    print('Performing Tukey Test on Reagents...')
                    tukey_reagents = pairwise_tukeyhsd(df['z_score'], df['Reagent_1_Short_Hand'])

                    # Forming the tukey dataframe.
                    tukey_df = pd.DataFrame(data=tukey_reagents._results_table.data[1:], columns=tukey_reagents._results_table.data[0])
                    tukey_df = tukey_df.loc[tukey_df['p-adj'] < 0.05].reset_index(drop=True)
                    if len(tukey_df) > 0:

                        # Finding top 15 outlier reagents.
                        unique_factors = np.unique(tukey_df['group1'].unique().tolist() + tukey_df['group2'].unique().tolist())
                        counts = [find_factor_count(tukey_df, x) for x in unique_factors]
                        top_15_factor_performance = pd.DataFrame({'factor' : unique_factors, 'count' : counts})
                        top_15_factor_performance = top_15_factor_performance.sort_values(by=['count'])[::-1][0:15]

                        # Finding the average z-score for each top outlier catalyst.
                        average_z_scores = [find_average_zscore_factor(df, 'Reagent_1_Short_Hand', x) for x in top_15_factor_performance['factor']]
                        top_15_factor_performance['average_z_score'] = average_z_scores
                        top_15_factor_performance = top_15_factor_performance[['factor', 'average_z_score']]
                        # Saving the dataframe.
                        if time_data is not None:
                            name = reaction_type + '_homo_hydrogenation_top_15_outlier_reagents_' + str(time_data) + '.csv'
                        elif remove_negative_rxns == True:
                            name = reaction_type + '_homo_hydrogenation_top_15_outlier_catalysts_no_0%_yield.csv'
                        else:
                            name = reaction_type + '_homo_hydrogenation_top_15_outlier_reagents.csv'
                        top_15_factor_performance.to_csv(os.path.join(save_dir, name))
                        print(f'Saved outlier reagents dataframe to {os.path.join(save_dir, name)}.')
                    else:

                        print('No reagent was statistically better than any other.')
        else:
            print('Neither catalysts nor reagents were statistically significant.')
    
    else:
        print(f'Too few reactions in dataset. This dataset has {len(df)} reactions. We recommend having over 100 reactions.')
if __name__ == '__main__':
    main()
