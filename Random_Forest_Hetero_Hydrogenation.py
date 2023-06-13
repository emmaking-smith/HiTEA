'''
Running Random Forest for variable importances on heterogeneous hydrogenation data.
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
    parser.add_argument('--data_path', '-d', type=str, default='data/cleaned_datasets/hetero_hydrogenation.csv')
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

hetero_rxn_type_dict = {'alkene' : ['c12c(ncnc1N1CCN(CC1)C(OC(C)(C)C)=O)C=C(CC2)c1c(ccc2c1cnn2[C@@H]1OCCCC1)C',
                                   'N1(C(=N[C@](C(S1(=O)=O)=C)(C)c1scc(n1)NC(c1ccc(cn1)OC(F)F)=O)NC(OC(C)(C)C)=O)C',
                                   'C=1(C(CCCC1NC(C)=O)=O)C'],
                        
                        'deprotection' : ['c1c(cc2c(c1)C(=C([C@H](O2)c1ccc(cc1)OCc1ccccc1)c1ccc(cc1)OC)C)OC',
                                         'C1(CC(C1)OCc1ccccc1)N1CCOCC1',
                                         'N1(CC(CC1)(C(=O)OCC)C(=O)OCC)Cc1ccccc1',
                                         'n1c(nc2c(c1c1cnc(nc1)N)CCN2[C@]1(CCN(C1)C(=O)OCc1ccccc1)C)N1CCOCC1',
                                         'c1cc(c2c(c1OCc1ccccc1)c(n(cc2)C(OC(C)(C)C)=O)=O)F',
                                         'C1C2(CC1C2)N(Cc1ccccc1)Cc1ccccc1',
                                         '[nH]1c(nc2c(c1=O)C[C@H](CN2)CCc1sc(cc1C)C(=O)N[C@H](C(OCc1ccccc1)=O)CCc1nnnn1Cc1ccccc1)N',
                                         '[nH]1c(nc2c(c1=O)c(c[nH]2)CCc1ccc(cc1)C(=O)N[C@H](C(=O)OCc1ccccc1)CCc1[nH]nnn1)N',
                                         'n1c(nc2c(c1c1cnc(nc1)N)CCN2[C@]1(CCN(C1)Cc1ccccc1)C)N1CCOCC1'],
                        
                        'dearomatization' : ['c1c(cncc1O)C', 
                                             'c1(cc(ccc1)CO)CC(O)=O',
                                             'n1(nc(cc1)c1nccc(c1)C)C',
                                            'c1c(ccc(c1)C[C@H](C(OC)=O)NC(OC(C)(C)C)=O)O',
                                            'c12nc([nH]c2cccc1C(=O)O)C'],
                        
                        'other_FG' : ['C1([C@@H]([C@@H]([C@H](O1)COC(=O)c1ccccc1)OC(c1ccccc1)=O)OC(c1ccccc1)=O)N=[N+]=[N-]',
                                     'c1(sc2c(c1)c(nn2C(C)OCC)Br)C(=O)OCC',
                                     'c1(cc(c2c(c1)c([nH]cc2)=O)[N+](=O)[O-])c1c(nnn1C)C',
                                     'c1cc2n(c(c1)C)n(c(c2C#N)=O)Cc1ccccc1'],
                        
                        'alkyne' : ['c12c(ncc(c1)C#Cc1sc(cc1)C(=O)OCC)nc([nH]c2=O)NC(C(C)(C)C)=O']
                        }

def filter_out_rxn_types(df, rxn_type):
    rxns = hetero_rxn_type_dict[rxn_type]
    filtered_df = pd.DataFrame()
    for i in range(len(df)):
        rxn = df.loc[i, 'Reactant_1_SMILES']
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
    df = df.drop(columns=['Reactant_1_SMILES'])
    return df
    
def standard_analysis(df, save_dir, remove_negative_rxns=False):
    if remove_negative_rxns == True:
        df = df.loc[df['Product_Yield_PCT_Area_UV'] > 0]
    df = df[['Product_Yield_PCT_Area_UV', 'catalyst', 'Reactant_1_SMILES', 'ligand', 'Reagent_1_Short_Hand', 'Solvent_1_Name', 'Reaction_T']]
    df = df.reset_index(drop=True)
    df = cat_lig_column(df)
    for k in list(hetero_rxn_type_dict.keys()):
        df_i = filter_out_rxn_types(df, k)
        print(f'{k} has {len(df_i)} rxns.')
        if len(df_i) > 100:
            df_setup = workflow(df_i, k)
            r = random_forest(df_setup)
            if len(r) > 0:
                if remove_negative_rxns == True:
                    name = 'Hetero Hydrog No 0% Yield Reactions : ' + k
                    save_name = 'hetero_hydrog_no_0%_yield_reactions_' + k + '.csv'
                else:
                    name = 'Hetero Hydrog : ' + k
                    save_name = 'hetero_hydrog_' + k + '.csv'
                rf_imp_plot(r, feat_dict, name)
                
                r.to_csv(os.path.join(save_dir, save_name))
            else:
                print('Random forest variable importance analysis failed to find any variables of high importance.')
        else:
            print(f'Too few reactions in {k} dataset. We recommend having over 100 reactions.')
            
    return 0
    
    
def temporal_changes(df, year, save_dir):
    df = df[['Year', 'Product_Yield_PCT_Area_UV', 'catalyst', 'Reactant_1_SMILES', 'ligand', 'Reagent_1_Short_Hand', 'Solvent_1_Name', 'Reaction_T']]
    df = df.reset_index(drop=True)
    df = cat_lig_column(df)
    df = df.loc[df['Year'] == year].reset_index(drop=True)
    for k in list(hetero_rxn_type_dict.keys()):
        df_i = filter_out_rxn_types(df, k)
        print(f'{year}, {k} has {len(df_i)} reactions.')
        if len(df_i) > 100:
            df_setup = workflow(df_i, k)
            r = random_forest(df_setup)
            if len(r) > 0:
                rf_imp_plot(r, feat_dict, 'Hetero Hydrog ' + str(year) + ' Data Only: ' + k)
                r.to_csv(os.path.join(save_dir, 'hetero_hydrog_' + k + '_' + str(year) + '.csv'))
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
    df.loc[pd.isna(df['Reagent_1_Short_Hand']) == True, ['Reagent_1_Short_Hand']] = 'None'
    df.loc[df['catalyst'] == 'RaNi 4100', ['catalyst']] = '[Ni]'
    df.loc[df['catalyst'] == 'RaNi 4200', ['catalyst']] = '[Ni]'
    df = df.loc[pd.isna(df['catalyst']) == False]

   
    if time_data is not None:
        temporal_changes(df, time_data, save_dir)
    else:
        standard_analysis(df, save_dir, remove_negative_rxns=remove_negative_rxns)

if __name__ == '__main__':
    main()
        
        
    
        
    
    
    
    
    
