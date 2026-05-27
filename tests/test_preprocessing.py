import pandas as pd
import pytest
from src.data_loading import parse_geo,load_annotation
from src.preprocessing import select_samples,gene_matrix

def test_real_source_integrity():
    matrix,meta=parse_geo(__import__('pathlib').Path('data/raw/GSE45827_series_matrix.txt.gz'))
    assert matrix.shape==(29873,155)
    selected,eligible=select_samples(matrix,meta)
    assert selected.shape==(29873,59)
    assert eligible.groupby(['subtype','batch']).size().groupby(level=0).sum().to_dict()=={'Luminal A':29,'Luminal B':30}
    genes,mapping=gene_matrix(selected,load_annotation(__import__('pathlib').Path('data/raw/GPL570.annot.gz')))
    assert genes.shape==(13823,59)
    assert not genes.index.duplicated().any()
    assert {'ESR1','PGR','MKI67'}.issubset(genes.index)

def test_duplicate_sample_fails():
    x=pd.DataFrame({'S1':[5.], 'S2':[6.]},index=['probe'])
    m=pd.DataFrame({'sample_id':['S1','S1'],'subtype':['Luminal A','Luminal B'],'source':['Tumor Sample']*2})
    with pytest.raises(ValueError,match='Duplicate'):
        select_samples(x,m)
