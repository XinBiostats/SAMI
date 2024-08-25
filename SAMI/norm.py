import numpy as np
import pandas as pd
import re
import os

def SumNorm(x, c=1): 
    if np.any(np.nansum(x, axis=1, keepdims=True)==0):
        print('Warning: The denominator array contains 0')
        
    out = c*x/np.nansum(x, axis=1, keepdims=True)
    
    return out

def MedianNorm(x): 
    if np.any(np.nanmedian(x, axis=1, keepdims=True)==0):
        print('Warning: The denominator array contains 0')
        
    out = x/np.nanmedian(x, axis=1, keepdims=True)
    
    return out

def CompNorm(x, ref_idx, c=1):
 
    ref_col = x[:,ref_idx][:, np.newaxis]
    
    if np.any(ref_col==0):
        print('Warning: The denominator array contains 0')
        
    out = c*x/ref_col
    
    return out

def SamplePQN(x, ref_smpl):
    
    if np.any(np.nanmedian(x/ref_smpl,axis=1,keepdims=True)==0):
        print('Warning: The denominator array contains 0')
        
    out = x/np.nanmedian(x/ref_smpl,axis=1,keepdims=True)
    
    return out

def GroupPQN(x,ref_smpl):

    return SamplePQN(x,np.nanmean(ref_smpl,axis=0))

def SpecNorm(x,norm_vec=None):

    if norm_vec is None:
        norm_vec = np.ones(x.shape[0])
        print("No sample specific information were given, all set to 1.0")
    else:
        norm_vec = np.array(norm_vec)
        
    if len(norm_vec.shape)==2: # to make sure norm_vec is a column vector
        out = x/norm_vec
    else:
        out = x/norm_vec[:,np.newaxis]
        
    return out

def GroupMedianPQN(x):

    grouped_medians_list = []
    for group in unique_groups:
        group_indices = np.where(np.array(group_labels) == group)[0]

        group_data = x[group_indices]

        median_values = np.median(group_data, axis=0)

        grouped_medians_list.append(median_values)

    grouped_medians = np.vstack(grouped_medians_list)

    overall_median = np.nanmedian(grouped_medians,axis=0)
    
    if np.any(grouped_medians==0):
        print('Warning: The denominator array contains 0')
    factor = overall_median/grouped_medians

    factor_dict = {}
    for i, row in enumerate(factor):
        factor_dict[unique_groups[i]] = row

    out = x.copy().astype(float)
    for group in unique_groups:
        group_indices = np.where(np.array(group_labels) == group)[0]
        out[group_indices] = x[group_indices]*factor_dict[group]

    return out

def LogTrans1(x, log_base=2):
    min_val = abs(x[x.nonzero()]).min()/100
    if log_base == 10:
        out = np.log10((x + np.sqrt(x**2 + min_val**2))/2)
    elif log_base == 2:
        out = np.log2((x + np.sqrt(x**2 + min_val**2))/2)
    
    return out 

def LogTrans2(x, log_base=2):
    
    if log_base == 10:
        out = np.log10(x + 1)
    elif log_base == 2:
        out = np.log2(x + 1)
    return out 

def SquareRootTrans(x):
    min_val = abs(x[x.nonzero()]).min()/10
    return np.sqrt((x + np.sqrt(x**2 + min_val**2))/2)

def CubeRootTrans(x):
    min_val = abs(x[x.nonzero()]).min()/10
    return ((x + np.sqrt(x**2 + min_val**2))/2)**(1/3)

def AutoNorm(x):
    return (x-np.nanmean(x,axis=0))/np.nanstd(x,axis=0)

def RangeNorm(x):
    return (x-np.nanmean(x,axis=0))/(np.nanmax(x,axis=0)-np.nanmin(x,axis=0))

def MeanCenter(x):
    return x-np.nanmean(x,axis=0)


def Normalization(df, first_compound_idx, rowNorm=None, transNorm=None, scaleNorm=None, ref_sample_idx=None, ref_compound=None, norm_vec=None, group_column=None, c=1, log_base=2):
    #first_compound_idx = df.columns.get_loc(first_compound)
    df = df[df.iloc[:,first_compound_idx:].sum(axis=1)!=0].reset_index(drop=True)
    data = df.iloc[:,first_compound_idx:].values.astype(float)
    
    ### sampel-wise normalization
    if rowNorm is not None and rowNorm not in ['SumNorm', 'MedianNorm', 'CompNorm', 'SamplePQN', 'GroupPQN', 'SpecNorm', 'GroupMedianPQN']:
        raise ValueError("'rowNorm' should be None or one of those methods: 'SumNorm', 'MedianNorm', 'CompNorm', 'SamplePQN', 'GroupPQN', 'SpecNorm', 'GroupMedianPQN'.")
    
    elif rowNorm == 'SumNorm':
        data = SumNorm(data,c)
        print(f"Normalization to constant sum {c}.")
        
    elif rowNorm == 'MedianNorm':
        data = MedianNorm(data)
        print("Normalization to sample median.")
    
    elif rowNorm == 'CompNorm':
        if ref_compound is None:
            raise ValueError('Please provide reference compound.')
        else:
            ref_idx = df.columns.get_loc(ref_compound)-first_compound_idx
            data = CompNorm(data,ref_idx,c)
            print("Normalization by a reference feature.")
            
    elif rowNorm == 'SamplePQN':
        if ref_sample_idx is None:
            raise ValueError('Please provide reference sample index.')
        else:
            ref_smpl = data[ref_sample_idx,:]
            data = SamplePQN(data,ref_smpl)
            print("Probabilistic Quotient Normalization by a reference sample.")
            
    elif rowNorm == 'GroupPQN':
        if ref_sample_idx is None:
            raise ValueError('Please provide reference group samples index.')
        else:
            ref_smpl = data[ref_sample_idx,:]
            data = GroupPQN(data,ref_smpl)
            print("Probabilistic Quotient Normalization by a reference group.")
            
    elif rowNorm == 'SpecNorm':
        data = SpecNorm(data,norm_vec)
        print("Normalization by sample-specific factor.")
        
    elif rowNorm == 'GroupMedianPQN':
        if group_column is None:
            raise ValueError('Please provide the name of the column indicating different groups.')
        else:
            group_labels = df[group_column].tolist()
            unique_groups = np.unique(group_labels)
            data = GroupMedianPQN(data)
            print("Probabilistic Quotient Normalization by a customized sample based on median of each group.")
           
    else:
        data = data
        print("Normalization: N/A.")
        
    # if the reference by feature, the feature column should be removed, since it is all 1
    if rowNorm == 'CompNorm' and ref_compound is not None:
        ref_idx = df.columns.get_loc(ref_compound)-first_compound_idx
        data = np.delete(data, ref_idx, axis=1)
    
    ### Transformation
    if transNorm is not None and transNorm not in ['LogTrans1', 'LogTrans2', 'SquareRootTrans', 'CubeRootTrans']:
        raise ValueError("'transNorm' should be None or one of those methods: 'LogTrans', 'SquareRootTrans', 'CubeRootTrans'.")
        
    elif transNorm == 'LogTrans1':
        data = LogTrans1(data, log_base)
        if log_base == 10:
            print("Log10 Transformation.")
        elif log_base == 2:
            print("Log2 Transformation.")
            
    elif transNorm == 'LogTrans2':
        data = LogTrans2(data, log_base)
        if log_base == 10:
            print("Log10 Transformation.")
        elif log_base == 2:
            print("Log2 Transformation.")
        
    elif transNorm == 'SquareRootTrans':
        data = SquareRootTrans(data)
        print("Square Root Transformation.")
        
    elif transNorm == 'CubeRootTrans':
        data = CubeRootTrans(data)
        print("Cubic Root Transformation.")
    
    else:
        data = data
        print("Transformation: N/A.")
        
    ### Scaling
    if scaleNorm is not None and scaleNorm not in ['MeanCenter', 'AutoNorm', 'ParetoNorm', 'RangeNorm']:
        raise ValueError("'scaleNorm' should be None or one of those methods: 'MeanCenter', 'AutoNorm', 'ParetoNorm', 'RangeNorm'.")
        
    elif scaleNorm == 'MeanCenter':
        data = MeanCenter(data)
        print("Mean Centering.")
    
    elif scaleNorm == 'AutoNorm':
        data = AutoNorm(data)
        print("Autoscaling.")
    
    elif scaleNorm == 'ParetoNorm':
        data = ParetoNorm(data)
        print("Pareto Scaling.")
    
    elif scaleNorm == 'RangeNorm':
        data = RangeNorm(data)
        print("Range Scaling.")
        
    else:
        data = data
        print("Scaling: N/A.")
    
    new_df = df.copy()
    
    if rowNorm == 'CompNorm' and ref_compound is not None:
        new_df = new_df.drop(columns = ref_compound)
        
    new_df.iloc[:,first_compound_idx:] = data
    
    return new_df

def create_norm_dataset(data_path,pattern):
    for file in os.listdir(data_path):
        name, ext = os.path.splitext(file)
        if re.match(pattern,file):
            print('-----------------------------')
            df = pd.read_csv(os.path.join(data_path,f"{file}"))
            df_norm = Normalization(df, first_compound_idx=3, rowNorm='SumNorm', transNorm='LogTrans2', c=1, log_base=2)
            name = f"{name}_norm"
            df_norm.to_csv(os.path.join(data_path,f"{name}{ext}"),index=False)
            print(f'{name}{ext} is created.')
    