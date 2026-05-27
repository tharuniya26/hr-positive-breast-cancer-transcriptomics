import numpy as np
import pandas as pd
import pytest
from src.quality_control import qc_metrics,pca

def test_missing_data_fails():
    x=pd.DataFrame([[1.,np.nan]],index=['g'],columns=['a','b'])
    m=pd.DataFrame({'sample_id':['a','b'],'subtype':['Luminal A','Luminal B'],'batch':['1','1']})
    with pytest.raises(ValueError,match='Missing'):
        qc_metrics(x,m)

def test_pca_orientation_and_dimensions():
    x=pd.DataFrame([[1,2,3,4],[4,3,2,1],[1,3,1,3]],index=['g1','g2','g3'],columns=list('abcd'))
    m=pd.DataFrame({'sample_id':list('abcd'),'subtype':['Luminal A']*2+['Luminal B']*2,'batch':['1']*4})
    a,ev=pca(x,m);b,_=pca(x,m)
    assert a.shape==(4,5)
    assert np.allclose(a[['PC1','PC2']],b[['PC1','PC2']])
    assert np.all((ev>=0)&(ev<=1))
