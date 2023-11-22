import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt



## dataset managing

# splitting of dataset int H and U sets according to "ishealthy" mask
def split_dataset(dataset,ishealthy):
    nonzero_rows = (dataset.sum(axis=1)>0)  
    Hds = dataset.loc[nonzero_rows].loc[:, ishealthy] # healthy view
    Uds = dataset.loc[nonzero_rows].loc[:,~ishealthy] # unhealthy view
    return Hds,Uds

# split dataset in balanced and unbalanced datasets optionally apply threshold
def split_dataset2(dataset,ishealthy,threshold=None):

    # apply threshold on functions relevance
    if threshold is not None:
        # apply threshold and retain only surviving species
        dataset = dataset[dataset>threshold].fillna(0)

    # remove unused functions
    #dataset = dataset.loc[dataset.sum(axis=1)>0]

    # splitting of dataset  
    Hds      = dataset.loc[:, ishealthy] # healthy view
    Uds_ubal = dataset.loc[:,~ishealthy] # unhealthy view

    # selcetion of balanced subset of unhealthy samples
    NH = ( ishealthy).sum()  # number of   healty samples
    NU = (~ishealthy).sum()  # number of unhealty samples
    samples_selection = np.random.choice(NU,NH,replace=False)
    samples_selection.sort()
    Uds_bal = Uds_ubal.iloc[:,samples_selection]

    return Hds, Uds_bal, Uds_ubal

# balanced randomized splitting
def split_dataset_bal(dataset,ishealthy,threshold=None):
    Hds, Uds_bal, Uds_ubal = split_dataset2(dataset,ishealthy)
    return Hds, Uds_bal

# compue threshold value with standard procedure
def compute_threshold(dataset, Nmesh=150, return_functions=False,logder=True):

    # instantiate vectorized "average alpha div" function
    average_alpha_div = np.vectorize(lambda tr:(dataset > tr).sum(axis=0).mean() )

    # define domain extrema
    vmin = dataset.min().min()+1e-9
    vmax = dataset.max().max()

    # compute function values
    threshold = np.logspace(np.log10(vmin),np.log10(vmax),Nmesh)
    alpha = average_alpha_div(threshold)

    # compute derivative values
    if logder:
        threshold2    = 10**((np.log10(threshold[1:])+np.log10(threshold[:-1]))/2)
        dlogthreshold =      (np.log10(threshold[1:])-np.log10(threshold[:-1]))
        dalpha        =      (             alpha[1:] -             alpha[:-1] )
        diff          = dalpha/dlogthreshold
    else:
        threshold2 = (threshold[1:]+threshold[:-1])/2
        dthreshold = (threshold[1:]-threshold[:-1])
        dalpha     = (    alpha[1:]-    alpha[:-1])
        diff = dalpha/dthreshold

    # smooth out diff
    diff = pd.Series(diff).rolling(int(np.ceil(0.05*Nmesh)),center=True).mean().fillna(pd.Series(diff)).values

    # compute threshold value
    thr = threshold2[np.argmin(diff)]

    # returns
    if return_functions:
        return thr, (threshold,alpha), (threshold2,diff)
    else:
        return thr




## bins managing

# logbinning within global min & max
def make_bins(Hds,Uds,nbins=30):

    # define domain extrema
    vmin = min(Hds[Hds>0].min().min(),Uds[Uds>0].min().min())
    vmax = max(Hds[Hds>0].max().max(),Hds[Hds>0].max().max())

    # define bins position
    bins = np.logspace(np.log10(vmin),np.log10(vmax),nbins+1)

    return bins

# reconstruct actual bin edges from bins centers
def reconstruct_bins(bins_center, logbins=True):
    if logbins:
        bins_exp  = np.log10(bins_center)
        delta_exp = bins_exp[1]-bins_exp[0]
        bins      = 10**np.concatenate([ [bins_exp[0]-delta_exp/2], bins_exp+delta_exp/2 ])
    else:
        delta      = bins_center[1]-bins_center[0]
        bins      = np.concatenate([ [bins_center[0]-delta/2], bins_center+delta/2 ])

    return bins





## normalizations

# normalizes as a probability distribution
def normalize_vec(x):
    s=sum(x)
    if s==0: return x*0
    else   : return x/s


# rescale the values so that the distributions along axis=axis are 0-mean 1-variance (in logvalues or values)
def rescale_dataset(ds,axis=0,log_scaling=True):
    if log_scaling:
        logds = np.log10(ds[ds>0])
        return 10**(  (logds - np.expand_dims( logds.median(axis=axis).values, axis=axis ) )/np.expand_dims( logds.std(axis=axis).values, axis=axis )  )
    else:
        return (ds - np.expand_dims( ds.median(axis=axis).values, axis=axis ) )/np.expand_dims( ds.std(axis=axis).values, axis=axis )





## SVD projection tools

def compute_svd_decomposition(ds,standardize=False):
    m = ds.mean(axis=1)
    s = ds.std(axis=1)
    if standardize:
        ds = ((ds - m.values[:,None])/s.values[:,None])
        ds.loc[s==0] = 0
    else:
        ds = (ds - m.values[:,None])

    U, S, Vh = np.linalg.svd(ds, full_matrices=False)
    SVD = {"m":m, "U":pd.DataFrame(U,index=ds.index), "S":S ,"Vh":pd.DataFrame(Vh,index=ds.columns),"s":s}

    return SVD

def project_on_subspace(ds,SVD,n,svd_basis=True):
    Ured = SVD["U"].iloc[:,:n]  # basis truncated to n elements
    if svd_basis:
        ds_red = Ured.T @ (ds - SVD["m"].values[:,None])
    else:
        ds_red = Ured @ (Ured.T @ (ds - SVD["m"].values[:,None]))

    return ds_red

def reconstruct_from_subspace(ds_red,SVD,add_mean=True):
    n = len(ds_red)             # dimension of subspace
    Ured = SVD["U"].iloc[:,:n]  # basis truncated to n elements

    if   len(ds_red.shape)==1: ds_reconstructed = (SVD["m"]                if add_mean else 0) + Ured @ ds_red
    elif len(ds_red.shape)==2: ds_reconstructed = (SVD["m"].values[:,None] if add_mean else 0) + Ured @ ds_red
        
    return ds_reconstructed







## other functions

# computes asymmetric standard deviation (geometric average of only positive/negative variations)
def asymmetric_std(data,axis=1):

    # contracts along array dimension
    if len(data.shape)==1:
        m = data.mean()
        delta = data - m
        n_pos = (delta>0).sum() + 0.5*(delta==0).sum()
        n_neg = len(data)-n_pos
        std_neg = np.sqrt( (np.minimum(0,delta)**2).sum()/np.maximum(1,n_neg) )
        std_pos = np.sqrt( (np.maximum(0,delta)**2).sum()/np.maximum(1,n_pos) )

    # contracts along axis=axis of dataframe
    else:
        m = data.mean(axis=axis)
        delta = data - m.values[:,None]
        n_pos = (delta>0).sum(axis=axis) + 0.5*(delta==0).sum(axis=axis)
        n_neg = delta.shape[1]-n_pos
        std_neg = np.sqrt( (np.minimum(0,delta)**2).sum(axis=axis)/np.maximum(1,n_neg) )
        std_pos = np.sqrt( (np.maximum(0,delta)**2).sum(axis=axis)/np.maximum(1,n_pos) )
    
    return std_neg,std_pos

# computes rolling average of each column (window over rows enumeration)
def rolling_average(dexpr,window=10):
    win_type = "blackmanharris"
    effective_amplitude = 0.277                     # 2*std of blackman-harris window (euristic)
    eff_window = round(window/effective_amplitude)  # rescale amplitude so that there is a corresoponding between linear binning and blackman-harris binning
    return dexpr.rolling(window=eff_window,center=True,win_type=win_type,closed="both",min_periods=1).mean()

# select only elements of upper triangular matrix (no diag)
def triu_elements(A):
    L = len(A)
    triu_mask = np.triu(np.full(shape=(L,L),fill_value=1),k=1).astype(bool)
    try: return A.values[triu_mask]
    except: return A[triu_mask]

# weight
def compute_weights(ds):
    mad = ds.sum(axis=1)
    w   = triu_elements(mad.values[None,:] * mad.values[:,None])
    return w


# binwidth indicator (if ax is logscaled binw is expegted to be logwidth of bin)
def add_binw_patch(ax,binw,axes_x0=0.05,axes_y0=0.03,axes_yw=0.03,axes_xw=0.02,axis=0,label=True): # defaut value sin axis coordinate system

    def make_binw_path(axis=0):
        vertices = np.array([(0,0),(1,0),(1,-0.5),(1,0.5),(0,-0.5),(0,0.5)])
        codes    = [1    ,2       ,1        ,2       ,1     ,2    ]  # 1=MOVETO, 2=LINETO
        rotation = mpl.transforms.Affine2D().rotate_around(x=0, y=0, theta=np.pi/2*axis) # identity if axis=0, else pi/2 rotation
        return mpl.path.Path(vertices=vertices, codes=codes,closed=False, readonly=False).transformed(rotation)

    def compute_x0_position(ax,axes_x0): # in scaled data coordinates
        xlim = ax.get_xlim()                              # triggers imit computtions
        tr_xlim = ax.transScale.transform([[xlim[0],0],[xlim[1],0]])[:,0] # xlim in scaled coordinates
        tr_xspan = tr_xlim[1]-tr_xlim[0]                  # width of x span
        x0 = tr_xlim[0] + tr_xspan*axes_x0                # where to position the origin of the bin indicator
        return x0
    
    def compute_y0_position(ax,axes_y0): # in scaled data coordinates
        ylim = ax.get_ylim()                              # triggers imit computtions
        tr_ylim = ax.transScale.transform([[0,ylim[0]],[0,ylim[1]]])[:,1] # ylim in scaled coordinates
        tr_yspan = tr_ylim[1]-tr_ylim[0]                  # width of y span
        y0 = tr_ylim[0] + tr_yspan*axes_y0                # where to position the origin of the bin indicator
        return y0

    assert axis in [0,1]
    if axis==0:
        x0 = compute_x0_position(ax,axes_x0=axes_x0) # compute x0 in scaled data coordinates (the xlims need to be fixed after this point)

        # coompute transformation that sends the "unit" binw patch in to the proper patch (x drawn on data layer, y drawn on axes layer)
        x_transform = mpl.transforms.Affine2D().scale(sx=binw,sy=0      ).translate(tx=x0,ty=0      ) + ax.transLimits + ax.transAxes
        y_transform = mpl.transforms.Affine2D().scale(sx=0   ,sy=axes_yw).translate(tx=0 ,ty=axes_y0) + ax.transAxes
        place = mpl.transforms.blended_transform_factory(x_transform, y_transform)
    if axis==1:
        y0 = compute_y0_position(ax,axes_y0=axes_y0) # compute x0 in scaled data coordinates (the xlims need to be fixed after this point)

        # coompute transformation that sends the "unit" binw patch in to the proper patch (x drawn on data layer, y drawn on axes layer)
        x_transform = mpl.transforms.Affine2D().scale(sx=axes_xw,sy=0   ).translate(tx=axes_x0,ty=0 ) + ax.transAxes
        y_transform = mpl.transforms.Affine2D().scale(sx=0      ,sy=binw).translate(tx=0      ,ty=y0) + ax.transLimits + ax.transAxes
        place = mpl.transforms.blended_transform_factory(x_transform, y_transform)

    # draw the patch
    patch = mpl.patches.PathPatch(make_binw_path(axis=axis),lw=2,zorder=2.5,transform=place)
    patch_handler = ax.add_patch(patch)
    
    # set label for binwidth indicator
    if label:
        import warnings
        with warnings.catch_warnings(): 
            warnings.simplefilter("ignore")
            proxy_artist = ax.errorbar([np.nan],[np.nan],xerr=[1], capsize=5,linewidth=2.5,capthick=2.5,color="black",label="binning",linestyle="")
        return patch_handler, proxy_artist
    
    else:
        return patch_handler
    

# extract subset of functions corresponding to some groups
# each function is weighted by its degree of specificity (1/#number of groups connected to it)
def functions_subset(groups,fg_annot,specific=True):
    if type(groups) is str: groups = [groups]  # common interface for scalar anf list input.

    if specific:
        mask = np.any( np.isclose(fg_annot[groups],1),axis=1)
    else:
        mask = np.any(~np.isclose(fg_annot[groups],0),axis=1)

    return fg_annot[groups][mask].sum(axis=1) 

