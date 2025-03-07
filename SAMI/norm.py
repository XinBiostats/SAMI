import numpy as np
import pandas as pd
import re
import os

def SumNorm(x, c=1): 
    """
    Performs sum normalization by scaling each sample to a constant sum.

    Parameters
    ----------
    x : np.ndarray
        Input matrix where rows represent samples and columns represent features.
    c : float, optional
        The target sum for normalization. Default is 1.

    Returns
    -------
    np.ndarray
        Normalized matrix with row sums equal to `c`.
    """
    if np.any(np.nansum(x, axis=1, keepdims=True)==0):
        print('Warning: The denominator array contains 0')
        
    out = c*x/np.nansum(x, axis=1, keepdims=True)
    
    return out

def MedianNorm(x):
    """
    Normalizes each sample by dividing by its median.

    Parameters
    ----------
    x : np.ndarray
        Input matrix.

    Returns
    -------
    np.ndarray
        Normalized matrix where each sample's median is scaled to 1.
    """
    if np.any(np.nanmedian(x, axis=1, keepdims=True)==0):
        print('Warning: The denominator array contains 0')
        
    out = x/np.nanmedian(x, axis=1, keepdims=True)
    
    return out

def CompNorm(x, ref_idx, c=1):
    """
    Performs normalization using a reference feature.

    Parameters
    ----------
    x : np.ndarray
        Input matrix.
    ref_idx : int
        Index of the reference feature.
    c : float, optional
        Scaling factor. Default is 1.

    Returns
    -------
    np.ndarray
        Normalized matrix where each value is divided by the reference feature.
    """
    ref_col = x[:,ref_idx][:, np.newaxis]
    
    if np.any(ref_col==0):
        print('Warning: The denominator array contains 0')
        
    out = c*x/ref_col
    
    return out

def SamplePQN(x, ref_smpl):
    """
    Applies Probabilistic Quotient Normalization (PQN) using a reference sample.

    Parameters
    ----------
    x : np.ndarray
        Input matrix.
    ref_smpl : np.ndarray
        Reference sample.

    Returns
    -------
    np.ndarray
        PQN normalized matrix.
    """    
    if np.any(np.nanmedian(x/ref_smpl,axis=1,keepdims=True)==0):
        print('Warning: The denominator array contains 0')
        
    out = x/np.nanmedian(x/ref_smpl,axis=1,keepdims=True)
    
    return out

def GroupPQN(x,ref_smpl):
    """
    Performs PQN using the mean of a reference group.

    Parameters
    ----------
    x : np.ndarray
        Input matrix.
    ref_smpl : np.ndarray
        Reference group samples.

    Returns
    -------
    np.ndarray
        Group PQN normalized matrix.
    """
    return SamplePQN(x,np.nanmean(ref_smpl,axis=0))

def SpecNorm(x,norm_vec=None):
    """
    Performs sample-specific normalization using an external normalization vector.

    Parameters
    ----------
    x : np.ndarray
        Input matrix.
    norm_vec : np.ndarray, optional
        Vector containing normalization factors for each sample. If None, defaults to 1.

    Returns
    -------
    np.ndarray
        Normalized matrix.
    """
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
    """
    Performs Probabilistic Quotient Normalization (PQN) using group-specific medians.

    Parameters
    ----------
    x : np.ndarray
        Input matrix where rows represent samples and columns represent features.
    group_labels : list or np.ndarray
        List of group labels corresponding to each row in `x`.

    Returns
    -------
    np.ndarray
        PQN normalized matrix where each sample is normalized by the median of its group.

    Notes
    -----
    - This method calculates the median for each group and scales individual samples accordingly.
    - The median across all group medians is used as a reference.
    - If any group median contains zeros, a warning is printed.
    """
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
    """
    Applies log transformation with a shift to avoid undefined values for zero or negative numbers.

    Parameters
    ----------
    x : np.ndarray
        Input matrix.
    log_base : int, optional
        Logarithm base (2 or 10). Default is 2.

    Returns
    -------
    np.ndarray
        Log-transformed matrix.

    Notes
    -----
    - A small shift (min nonzero value / 100) is applied to prevent log(0) issues.
    - Uses the transformation: log_base((x + sqrt(x^2 + min_val^2)) / 2).
    """
    min_val = abs(x[x.nonzero()]).min()/100
    if log_base == 10:
        out = np.log10((x + np.sqrt(x**2 + min_val**2))/2)
    elif log_base == 2:
        out = np.log2((x + np.sqrt(x**2 + min_val**2))/2)
    
    return out 

def LogTrans2(x, log_base=2):
    """
    Applies log transformation with a pseudocount of 1.

    Parameters
    ----------
    x : np.ndarray
        Input matrix.
    log_base : int, optional
        Logarithm base (2 or 10). Default is 2.

    Returns
    -------
    np.ndarray
        Log-transformed matrix.

    Notes
    -----
    - Uses the transformation: log_base(x + 1).
    """    
    if log_base == 10:
        out = np.log10(x + 1)
    elif log_base == 2:
        out = np.log2(x + 1)
    return out 

def SquareRootTrans(x):
    """
    Applies square root transformation to stabilize variance.

    Parameters
    ----------
    x : np.ndarray
        Input matrix.

    Returns
    -------
    np.ndarray
        Square root-transformed matrix.

    Notes
    -----
    - A small shift (min nonzero value / 10) is applied to avoid instability at small values.
    """
    min_val = abs(x[x.nonzero()]).min()/10
    return np.sqrt((x + np.sqrt(x**2 + min_val**2))/2)

def CubeRootTrans(x):
    """
    Applies cube root transformation to reduce skewness.

    Parameters
    ----------
    x : np.ndarray
        Input matrix.

    Returns
    -------
    np.ndarray
        Cube root-transformed matrix.

    Notes
    -----
    - A small shift (min nonzero value / 10) is applied for numerical stability.
    """
    min_val = abs(x[x.nonzero()]).min()/10
    return ((x + np.sqrt(x**2 + min_val**2))/2)**(1/3)

def AutoNorm(x):
    """
    Performs autoscaling (z-score normalization).

    Parameters
    ----------
    x : np.ndarray
        Input matrix.

    Returns
    -------
    np.ndarray
        Z-score normalized matrix.

    Notes
    -----
    - The transformation: (x - mean) / standard deviation.
    - Each feature is scaled to have mean 0 and standard deviation 1.
    """
    return (x-np.nanmean(x,axis=0))/np.nanstd(x,axis=0)

def RangeNorm(x):
    """
    Normalizes data by centering and scaling to the range.

    Parameters
    ----------
    x : np.ndarray
        Input matrix.

    Returns
    -------
    np.ndarray
        Range-normalized matrix.

    Notes
    -----
    - Uses the transformation: (x - mean) / (max - min).
    """
    return (x-np.nanmean(x,axis=0))/(np.nanmax(x,axis=0)-np.nanmin(x,axis=0))

def MeanCenter(x):
    """
    Performs mean centering by subtracting the mean of each feature.

    Parameters
    ----------
    x : np.ndarray
        Input matrix.

    Returns
    -------
    np.ndarray
        Mean-centered matrix.

    Notes
    -----
    - The transformation: x - mean.
    - Each feature will have a mean of zero.
    """
    return x-np.nanmean(x,axis=0)


def Normalization(df, first_compound_idx, rowNorm=None, transNorm=None, scaleNorm=None, ref_sample_idx=None, ref_compound=None, norm_vec=None, group_column=None, c=1, log_base=2):
    """
    Performs data normalization, transformation, and scaling on a given dataset.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe containing metabolomics/lipidomics data.
    first_compound_idx : int
        Index of the first compound column in `df` (i.e., the first feature column).
    rowNorm : str, optional
        Sample-wise normalization method. Options:
        - 'SumNorm': Normalize to constant sum `c`.
        - 'MedianNorm': Normalize by median of each sample.
        - 'CompNorm': Normalize by a specific reference compound.
        - 'SamplePQN': Probabilistic Quotient Normalization using a reference sample.
        - 'GroupPQN': PQN using a reference group of samples.
        - 'SpecNorm': Normalize using a given normalization vector.
        - 'GroupMedianPQN': PQN using median per group.
    transNorm : str, optional
        Transformation method. Options:
        - 'LogTrans1': Log transformation with offset to avoid undefined values.
        - 'LogTrans2': Log transformation with pseudocount (+1).
        - 'SquareRootTrans': Square root transformation.
        - 'CubeRootTrans': Cube root transformation.
    scaleNorm : str, optional
        Scaling method. Options:
        - 'MeanCenter': Mean centering.
        - 'AutoNorm': Autoscaling (z-score normalization).
        - 'ParetoNorm': Pareto scaling.
        - 'RangeNorm': Range scaling.
    ref_sample_idx : int, optional
        Index of the reference sample for `SamplePQN` or `GroupPQN`.
    ref_compound : str, optional
        Name of the reference compound for `CompNorm`.
    norm_vec : np.ndarray, optional
        Normalization vector for `SpecNorm`.
    group_column : str, optional
        Column name indicating groups (required for `GroupMedianPQN`).
    c : int, optional
        Constant for `SumNorm`. Default is 1.
    log_base : int, optional
        Logarithm base (2 or 10). Default is 2.

    Returns
    -------
    pd.DataFrame
        Normalized, transformed, and scaled dataframe.

    Notes
    -----
    - The function follows a sequential pipeline: **Normalization → Transformation → Scaling**.
    - If `rowNorm` is `CompNorm`, the reference compound column is removed after normalization.
    - If `group_column` is provided, `GroupMedianPQN` applies PQN using the median per group.
    - Ensures that feature columns containing all zeros are removed before processing.
    """
    
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
    """
    Applies normalization to all datasets in a specified directory that match a given filename pattern.

    Parameters
    ----------
    data_path : str
        Path to the directory containing the datasets.
    pattern : str
        Regular expression pattern to match filenames.

    Returns
    -------
    None
        Saves the normalized datasets as CSV files in the same directory.

    Notes
    -----
    - Uses `Normalization()` with `SumNorm` row normalization and `LogTrans2` transformation.
    - Adds `_norm` to the filename of the processed files.
    """
    for file in os.listdir(data_path):
        name, ext = os.path.splitext(file)
        if re.match(pattern,file):
            print('-----------------------------')
            df = pd.read_csv(os.path.join(data_path,f"{file}"))
            df_norm = Normalization(df, first_compound_idx=3, rowNorm='SumNorm', transNorm='LogTrans2', c=1, log_base=2)
            name = f"{name}_norm"
            df_norm.to_csv(os.path.join(data_path,f"{name}{ext}"),index=False)
            print(f'{name}{ext} is created.')
    