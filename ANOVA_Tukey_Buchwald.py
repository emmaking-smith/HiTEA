'''
Run an ANOVA-Tukey test on Buchwald data.
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
    parser.add_argument('--data_path', '-d', type=str, default='data/cleaned_datasets/buchwald.csv')
    parser.add_argument('--save_dir', '-s', type=str, default='.')
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

def pairs(df):
    pairs = []
    for i in range(len(df)):
        nuc = df.loc[i, 'nuc']
        hal = df.loc[i, 'halide']
        pairs.append(hal + '+' + nuc)
    df['pair'] = pairs
    df = df.drop(columns=['halide', 'nuc'])
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

buchwald_rxn_type_dict = {'br+sec_amine' : ['Brc1ccc2ncccc2c1+O=C1NCc2c(OCCCCN3CCNCC3)cccc21',
                                           'Brc1ccccn1+Cn1cc(Nc2nc(O[C@H]3C[C@H](N(C)C(=O)OC(C)(C)C)C3)c3ccn(COCC[Si](C)(C)C)c3n2)cn1',
                                           'CCc1ccn2ccnc2c1Br+CC(C)(C)OC(=O)N1CCCNCC1'],
                          
                          'cl+pri_amine' : ['CN(C(=O)OC(C)(C)C)[C@H]1C[C@H](Oc2nc(Cl)nc3[nH]cc(-c4ccccn4)c23)C1+Cn1cc(N)cn1',
                                           'Fc1c(Cl)nc(N2CCOCC2)nc1Cl+CC(C)(C)OC(=O)N1CC[C@](N)(CO)C1',
                                           'Cc1cc(C)c(CN2CCc3ccc(OS(=O)(=O)C(F)(F)F)c(Cl)c3C2=O)c(OCc2ccccc2)n1+CCN.Cl',
                                           'Cn1ccnc1Cl+Nc1cccnc1',
                                           'Nc1ccc(-c2cc(Cl)nc(N3CCOCC3)n2)cn1+CN(C(=O)OC(C)(C)C)[C@H]1C[C@H](N)C1',
                                           'CC(C)(C)OC(=O)N1C[C@H](COc2nc(Cl)nc3[nH]ccc23)[C@@H](C(F)(F)F)C1+Nc1cnoc1'],
                          
                          'cl+sec_amine' : ['Cn1cc(Nc2nc(Cl)nc3[nH]cnc23)cn1+CC(C)(C)OC(=O)N1C[C@@H]2CNC[C@@H]2C1',
                                           'Cc1ncc(Cl)nc1Cl+C1COCCN1',
                                           'Cc1cc(C)c(CN2CCc3ccc(OS(=O)(=O)C(F)(F)F)c(Cl)c3C2=O)c(OCc2ccccc2)n1+CCNC1CCOCC1'],
                          
                          'br+pri_amine' : ['Brc1ccc(N2CCOCC2)nc1+CC(C)(C)OC(=O)N1CC[C@@H](N)C1', 
                                            'CC(C)(C)OC(=O)N1Cc2c(Br)nn(C(=O)OC(C)(C)C)c2C1(C)C+Cn1cc(N)cn1',
                                           'Cn1ccnc1Br+Nc1cccnc1',
                                           'CC(C)(C)OC(=O)N1CC[C@@H](N)C1'],
                          
                          'br+amide' : ['Brc1cnsc1+CC(C)(C)OC(N)=O',
                                        'Brc1ccc2ncccc2c1+CC(C)(C)OC(=O)NNC(=O)OC(C)(C)C'],
                          
                          'br+aro_N' : ['O=C1OCCn2c1ccc(Br)c2=O+Cc1c[nH]cn1',
                                       'Brc1ccc2ncccc2c1+CCOC(=O)c1c[nH]nc1O'],
                          
                          'cl+pri_alc' : ['Cn1cc(Nc2nc(Cl)c3c(Cl)c[nH]c3n2)cn1+CC(C)(C)OC(=O)N1CC[C@@H](CO)C1'],
                          
                          'i+aro_N' : ['Cc1sc2nc(NCc3ccccc3)cnc2c1I+Cc1ncnc2[nH]ccc12'],
                          
                          'i+pri_amine' : ['Cn1cc(I)cn1+CCOC(=O)n1nc(N)c2c1C(C)(C)N(C(=O)OC(C)(C)C)C2',
                                          'Ic1cccnc1+Cn1ccnc1N'],
                          
                          'i+amide' : ['Clc1ccc(I)cc1+Cc1nnc2n1-c1ccccc1NC(=O)C2'],
                          
                          'br+pri_alc' : ['Cc1cc(Br)n(C)n1+CC(C)(C)OC(=O)CO',
                                         'Nc1ncncc1Br+CC(C)(C)OC(=O)N1CCC(CO)CC1'],
                          
                          'i+pri_alc' : ['Nc1ncncc1I+CC(C)(C)OC(=O)N1CCC(CO)CC1',
                                        'CC1=C(c2ccc(OC3CCCCO3)cc2)C(c2ccc(I)cc2)Oc2ccc(OC3CCCCO3)cc21+C[C@@H](CO)N1CC[C@@H](C)C1'],
                          
                          'cl+aro_N' : ['COc1ncccc1-c1cnc(N)c(O[C@H](C)c2cc(F)cnc2Cl)c1+c1c[nH]nn1'],
                          
                          'i+sec_alc' : ['CCOc1ccccc1I+CC(C)(C)OC(=O)N1CCC[C@@H](O)C1'],
                          
                          'br+sec_alc' : ['CCOc1ccccc1Br+CC(C)(C)OC(=O)N1CCC[C@@H](O)C1']
                        
                         }


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
    time_data = args.time_data
    remove_negative_rxns = args.remove_negative_rxns
    
    Path(save_dir).mkdir(parents=True, exist_ok=True)
    
    # Import data.
    df = pd.read_csv(data_path)
    if remove_negative_rxns == True:
        df = df.loc[df['Product_Yield_PCT_Area_UV'] > 0]
        
    df = df.loc[pd.isna(df['Reaction_T']) == False]
    df = df.loc[pd.isna(df['Reaction_Time_hrs']) == False].reset_index(drop=True)
    df = df.loc[df['Reagent_1_Short_Hand'] != 'Special']
    df = df.loc[pd.isna(df['catalyst']) == False]

    df.loc[df['Reagent_1_Short_Hand'] == 'KOPnt', ['Reagent_1_Short_Hand']] = 'KOtPn' 
    df.loc[df['Reagent_1_Short_Hand'] == 'LiOBut 1M in Hexanes', ['Reagent_1_Short_Hand']] = 'LiOBut' 
    
    if time_data is not None:
        df = df.loc[df['Year'] == time_data]
        
    df = df[['Product_Yield_PCT_Area_UV', 'PRODUCT_STRUCTURE', 'Solvent_1_Name', 'Reaction_T', 'halide', 'nuc', 
             'catalyst', 'ligand', 'Reagent_1_Short_Hand']]
    df = df.reset_index(drop=True)
    
    ### REMOVE AFTER ALPHA   
    print(len(df))
    idxs = list(df.index)
    to_include = 0.1
    keep_idxs = np.random.choice(idxs, int(np.ceil(len(idxs) * to_include) ), replace=False)
    df = df.loc[keep_idxs].reset_index(drop=True)
    print(len(df))
    ### END OF REMOVE
    
    if len(df) > 100:

        # Processing the data.
        df = cat_lig_column(df)
        df = pairs(df)

        df.loc[df['cat_lig'] == "c1(cc(cc(c1c1c(ccc(c1P(C1CCCCC1)C1CCCCC1)OC)OC)C(C)C)C(C)C)C(C)C.c1(c(cccc1)[Pd]Cl)CCN+", ['cat_lig']] = "BrettPhos Palladacycle+"
        df.loc[df['cat_lig'] == "C1(C)(C)c2c(c(P(c3ccccc3)c3ccccc3)ccc2)Oc2c1cccc2P(c1ccccc1)c1ccccc1.c1ccc(c(c1)[Pd+])c1ccccc1NC.[O-]S(C)(=O)=O+", ['cat_lig']] = "Xantphos Pd G4+"

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
                            name = 'buchwald_top_15_outlier_catalysts_' + str(time_data) + '.csv'
                        
                        elif remove_negative_rxns == True:
                            name = 'buchwald_top_15_outlier_catalysts_no_0%_yield' + '.csv'
                        else:
                            name = 'buchwald_top_15_outlier_catalysts.csv'
                        #top_15_factor_performance.to_csv(os.path.join(save_dir, name))
                        top_15_factor_performance.to_csv(os.path.join(save_dir, 'buchwald_10%_outlier_catalysts_trial_3.csv'))
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
                            name = 'buchwald_top_15_outlier_reagents_' + str(time_data) + '.csv'
                       
                        elif remove_negative_rxns == True:
                            name = 'buchwald_top_15_outlier_reagents_no_0%_yield' + '.csv'
                        else:
                            name = 'buchwald_top_15_outlier_reagents.csv'
                        #top_15_factor_performance.to_csv(os.path.join(save_dir, name))
                        top_15_factor_performance.to_csv(os.path.join(save_dir, 'buchwald_10%_outlier_reagents_trial_3.csv'))
                        print(f'Saved outlier reagents dataframe to {os.path.join(save_dir, name)}.')
                    else:

                        print('No reagent was statistically better than any other.')
        else:
            print('Neither catalysts nor reagents were statistically significant.')
    
    else:
        print(f'Too few reactions in dataset. This dataset has {len(df)} reactions. We recommend having over 100 reactions.')
if __name__ == '__main__':
    main()
