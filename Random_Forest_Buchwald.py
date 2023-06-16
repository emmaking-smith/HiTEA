'''
Running Random Forest for variable importances on Buchwald data.
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
    parser.add_argument('--data_path', '-d', type=str, default='data/cleaned_datasets/buchwald.csv')
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

def filter_out_rxn_types(df, rxn_type):
    rxns = buchwald_rxn_type_dict[rxn_type]
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
    df.loc[df['cat_lig'] == "c1(cc(cc(c1c1c(ccc(c1P(C1CCCCC1)C1CCCCC1)OC)OC)C(C)C)C(C)C)C(C)C.c1(c(cccc1)[Pd]Cl)CCN+", ['cat_lig']] = "BrettPhos Palladacycle+"
    df.loc[df['cat_lig'] == "C1(C)(C)c2c(c(P(c3ccccc3)c3ccccc3)ccc2)Oc2c1cccc2P(c1ccccc1)c1ccccc1.c1ccc(c(c1)[Pd+])c1ccccc1NC.[O-]S(C)(=O)=O+", ['cat_lig']] = "Xantphos Pd G4+"
    df = pairs(df)
    for k in list(buchwald_rxn_type_dict.keys()):
        df_i = filter_out_rxn_types(df, k)
        print(f'{k} has {len(df_i)} rxns.')
        if len(df_i) > 100:
            df_setup = workflow(df_i, k)
            r = random_forest(df_setup)
            if len(r) > 0:
                if remove_negative_rxns == True:
                    name = 'Buchwald No 0% Yield Reactions : ' + k
                    save_name = 'buchwald_no_0%_yield_reactions_' + k + '.csv'
                else:
                    name = 'Buchwald : ' + k
                    save_name = 'buchwald_' + k + '.csv'
                rf_imp_plot(r, feat_dict, name)
                
                #r.to_csv(os.path.join(save_dir, save_name))
                r.to_csv(os.path.join(save_dir, 'buchwald_' + k + '_10%_dataset_rf_importances_trial_1.csv'))
            else:
                print('Random forest variable importance analysis failed to find any variables of high importance.')
        else:
            print(f'Too few reactions in {k} dataset. We recommend having over 100 reactions.')
            
    return 0
    
    
def temporal_changes(df, year, save_dir):
    df = df[['Year', 'Product_Yield_PCT_Area_UV', 'catalyst', 'halide', 'nuc', 'ligand', 'Reagent_1_Short_Hand', 'Solvent_1_Name', 'Reaction_T']]
    df = df.reset_index(drop=True)
    df = cat_lig_column(df)
    df.loc[df['cat_lig'] == 'c1(cc(cc(c1c1c(ccc(c1P(C1CCCCC1)C1CCCCC1)OC)OC)C(C)C)C(C)C)C(C)C.c1(c(cccc1)[Pd]Cl)CCN+', ['cat_lig']] == 'BrettPhos Palladacycle+'
    df.loc[df['cat_lig'] == 'C1(C)(C)c2c(c(P(c3ccccc3)c3ccccc3)ccc2)Oc2c1cccc2P(c1ccccc1)c1ccccc1.c1ccc(c(c1)[Pd+])c1ccccc1NC.[O-]S(C)(=O)=O+', ['cat_lig']] == 'Xantphos Pd G4+'
    df = pairs(df)
    df = df.loc[df['Year'] == year].reset_index(drop=True)
    for k in list(buchwald_rxn_type_dict.keys()):
        df_i = filter_out_rxn_types(df, k)
        print(f'{year}, {k} has {len(df_i)} reactions.')
        if len(df_i) > 100:
            df_setup = workflow(df_i, k)
            r = random_forest(df_setup)
            if len(r) > 0:
                rf_imp_plot(r, feat_dict, 'Buchwald ' + str(year) + ' Data Only: ' + k)
                r.to_csv(os.path.join(save_dir, 'buchwald_' + k + '_' + str(year) + '.csv'))
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
    df = df.loc[df['Reagent_1_Short_Hand'] != 'Special']
    df = df.loc[pd.isna(df['catalyst']) == False]

    df.loc[df['Reagent_1_Short_Hand'] == 'KOPnt', ['Reagent_1_Short_Hand']] = 'KOtPn' 
    df.loc[df['Reagent_1_Short_Hand'] == 'LiOBut 1M in Hexanes', ['Reagent_1_Short_Hand']] = 'LiOBut'

    if time_data is not None:
        temporal_changes(df, time_data, save_dir)
    else:
        standard_analysis(df, save_dir, remove_negative_rxns=remove_negative_rxns)

if __name__ == '__main__':
    main()
        
        
    
        
    
    
    
    
    
