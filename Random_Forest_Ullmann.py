'''
Running Random Forest for variable importances on Ullmann data.
'''

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from matplotlib.pyplot import figure
import os
import argparse
import pickle
from pathlib import Path

# Load the feat_dict
with open('feature_dict.pickle', 'rb') as f:
    feat_dict = pickle.load(f)

def init_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_path', '-d', type=str, default='data/cleaned_datasets/ullmann.csv')
    parser.add_argument('--save_dir', '-s', type=str, default='.')
    parser.add_argument('--time_data', '-t', type=int, default=None, help='Are you using time data or not and if so, what year are you investigating? Default is no temporal analysis.')
    parser.add_argument('--remove_negative_rxns', '-n', type=bool, default=False, help='Do you wish to remove the 0% yield reactions (T/F)? Default is False, no 0% yield reactions removed.')
    return parser.parse_args()

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

ullmann_rxn_type_dict = {'i+aro_N' : ['Clc1ccc(I)cc1+CCOC(=O)c1c[nH]nc1O', 
                                      'CCOc1cc(Cl)c(C#N)cc1I+c1cn[nH]c1',
                                     'Cc1sc2nc(NCc3ccccc3)cnc2c1I+Clc1ncnc2[nH]ccc12'],
                        'i+pri_alc' : ['Nc1ncncc1I+CC(C)(C)OC(=O)N1CCC(CO)CC1',
                                      'CC1=C(c2ccc(OC3CCCCO3)cc2)C(c2ccc(I)cc2)Oc2ccc(OC3CCCCO3)cc21+C[C@@H](CO)N1CC[C@@H](C)C1'],
                        
                        'i+sec_alc' : ['COc1ccc(C2=C(C)c3ccc(OC)cc3O[C@@H]2c2ccc(I)cc2)cc1+CCCN1CC(O)C1',
                                      'CCOc1ccccc1I+CC(C)(C)OC(=O)N1CCC[C@@H](O)C1'],
                         
                         'br+pri_alc' : ['Cc1cc(Br)nn1C+C=CCO',
                                        'Nc1ncncc1Br+CC(C)(C)OC(=O)N1CCC(CO)CC1']
                        }

def filter_out_rxn_types(df, rxn_type):
    rxns = ullmann_rxn_type_dict[rxn_type]
    filtered_df = pd.DataFrame()
    for i in range(len(df)):
        rxn = df.loc[i, 'pair']
        if rxn in rxns:
            filtered_df = pd.concat((filtered_df, df[i:i+1]))
    return filtered_df.reset_index(drop=True)

def random_forest(df, drop_pairs=True):
    '''
    Finding the variable importances.
    Args:
        df (dataframe): The dataframe of the one-hot encoded dataset.
        drop_pairs (bool): True (default) - the importance of the reacting molecules/pairs 
                            are not shown in the variable importances.
    Returns:
        rf_importances (dataframe): All importances that are greater than 0.01.
    '''
    
    rf = RandomForestRegressor()

    train, test = train_test_split(df, test_size=0.15, random_state=9)
    train_input = train.loc[:, train.columns != 'Product_Yield_PCT_Area_UV']
    train_labels = train['Product_Yield_PCT_Area_UV']
    test_input = test.loc[:, test.columns != 'Product_Yield_PCT_Area_UV']
    test_labels = test['Product_Yield_PCT_Area_UV']
    
    rf.fit(train_input, train_labels)
    pred = rf.predict(test_input)
    
    importances = rf.feature_importances_
    
    feature_names = train_input.columns
    rf_importances = pd.DataFrame({'importances' : importances, 'feat' : feature_names})
    rf_importances = rf_importances.loc[rf_importances['importances'] > 0.01]
    rf_importances = rf_importances.reset_index(drop=True)
    if drop_pairs == True:
        # Removing pairs:
        drop_idxs = []
        for i in range(len(rf_importances)):
            match = re.search('pair_', rf_importances.loc[i, 'feat'])
            if match is not None:
                drop_idxs.append(i)
        rf_importances = rf_importances.drop(drop_idxs)
        rf_importances = rf_importances.reset_index(drop=True)
    
    return rf_importances

def get_feats(rf_importances, feat_dict):
    '''
    Updating the feature importances in the feature dictionary (feat_dict)
    Args:
        rf_importances (dataframe): The dataframe from random_forest()
        feat_dict (dictionary): The growing dictionary for the short hand catalyst/ligand/reagent/reactant names.
    Returns:
        feats (list): The list of features that are within the variable importances dataframe.
    '''
    feats = []
    for f in rf_importances['feat']:
        try:
            feats.append(feat_dict[f])
        except:
            feats.append(f)
    return feats

def rf_imp_plot(rf_importances, feat_dict, title_name):
    '''
    Creating the plots to visualize the variable importances.
    Args:
        rf_importances (dataframe): The dataframe from random_forest()
        feat_dict (dictionary): The growing dictionary for the short hand catalyst/ligand/reagent/reactant names.
        title_name (string): The name of the chart.
    Returns:
        NULL
    '''
    feats = get_feats(rf_importances, feat_dict)

    figure(figsize=(3.5, 3.5), dpi=150)

    plt.bar(feats, rf_importances['importances'])
    plt.title(r'Top RF Importances:' 
                '\n'
                r'' + title_name)
    plt.xlabel('Feature')
    plt.ylabel('Random Forest Importance')
    plt.xticks(rotation=90, size=7)
    plt.grid()
    plt.show()
    return
    
def workflow(df, rxn_type):
    df = filter_out_rxn_types(df, rxn_type)
    df = one_hot_solvent(df)
    df = one_hot_reagent(df)
    df = one_hot_cat_lig(df)
    df = df.drop(columns=['pair'])
    return df
    
def standard_analysis(df, save_dir, remove_negative_rxns=False):
    if remove_negative_rxns == True:
        df = df.loc[df['Product_Yield_PCT_Area_UV'] > 0]
    df = df[['Product_Yield_PCT_Area_UV', 'catalyst', 'halide', 'nuc', 'ligand', 'Reagent_1_Short_Hand', 'Solvent_1_Name', 'Reaction_T']]
    df = df.reset_index(drop=True)
    df = cat_lig_column(df)
    df = pairs(df)
    for k in list(ullmann_rxn_type_dict.keys()):
        df_i = filter_out_rxn_types(df, k)
        print(f'{k} has {len(df_i)} rxns.')
        if len(df_i) > 100:
            df_setup = workflow(df_i, k)
            r = random_forest(df_setup)
            if len(r) > 0:
                if remove_negative_rxns == True:
                    name = 'Ullmann No 0% Yield Reactions : ' + k
                    save_name = 'ullmann_no_0%_yield_reactions_' + k + '.csv'
                else:
                    name = 'Ullmann : ' + k
                    save_name = 'ullmann_' + k + '.csv'
                rf_imp_plot(r, feat_dict, name)
                
                r.to_csv(os.path.join(save_dir, save_name))
            else:
                print('Random forest variable importance analysis failed to find any variables of high importance.')
        else:
            print(f'Too few reactions in {k} dataset. We recommend having over 100 reactions.')
            
    return 0
    
    
def temporal_changes(df, year, save_dir):
    df = df[['Year', 'Product_Yield_PCT_Area_UV', 'catalyst', 'halide', 'nuc', 'ligand', 'Reagent_1_Short_Hand', 'Solvent_1_Name', 'Reaction_T']]
    df = df.reset_index(drop=True)
    df = cat_lig_column(df)
    df = pairs(df)
    df = df.loc[df['Year'] == year].reset_index(drop=True)
    for k in list(ullmann_rxn_type_dict.keys()):
        df_i = filter_out_rxn_types(df, k)
        print(f'{year}, {k} has {len(df_i)} reactions.')
        if len(df_i) > 100:
            df_setup = workflow(df_i, k)
            r = random_forest(df_setup)
            if len(r) > 0:
                rf_imp_plot(r, feat_dict, 'Ullmann ' + str(year) + ' Data Only: ' + k)
                r.to_csv(os.path.join(save_dir, 'ullmann_' + k + '_' + str(year) + '.csv'))
            else:
                print('Random forest found no variables of high importance.')
        else:
            print('Too few reactions in dataset. We recommend having over 100 reactions.')
    return 0

def main():
    args = init_args()
    data_path = args.data_path
    save_dir = args.save_dir
    time_data = args.time_data
    remove_negative_rxns = args.remove_negative_rxns
    
    Path(save_dir).mkdir(parents=True, exist_ok=True)
    
    # Import data.
    df = pd.read_csv(data_path)
    df = df.loc[pd.isna(df['Reaction_T']) == False]
    df = df.loc[pd.isna(df['Reaction_Time_hrs']) == False].reset_index(drop=True)
    df = df.loc[pd.isna(df['catalyst']) == False]

   
    if time_data is not None:
        temporal_changes(df, time_data, save_dir)
    else:
        standard_analysis(df, save_dir, remove_negative_rxns=remove_negative_rxns)

if __name__ == '__main__':
    main()
        
        
    
        
    
    
    
    
    
