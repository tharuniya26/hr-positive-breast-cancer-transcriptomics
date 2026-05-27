from pathlib import Path
import numpy as np
import pandas as pd
from src.differential_expression import fit_moderated
from src.pathway_analysis import overrepresentation,read_reactome
from src.visualization import workflow

def test_moderated_model_and_bh():
    rng=np.random.default_rng(10)
    ids=[f'S{i}' for i in range(20)]
    metadata=pd.DataFrame({'sample_id':ids,'subtype':['Luminal A']*10+['Luminal B']*10,
                           'batch':['1','2']*10})
    values=rng.normal(6,.4,(100,20));values[0,10:]+=3
    x=pd.DataFrame(values,index=[f'G{i}' for i in range(100)],columns=ids)
    out=fit_moderated(x,metadata).set_index('gene_symbol')
    assert out.loc['G0','log2_fold_change']>2
    assert out.loc['G0','adjusted_p_value']<.05
    assert out.adjusted_p_value.ge(out.p_value).all()
    assert out.adjusted_p_value.between(0,1).all()

def test_ora_universe_and_figure(tmp_path):
    de=pd.DataFrame({'gene_symbol':[f'G{i}' for i in range(30)],
                     'adjusted_p_value':[.001]*10+[.9]*20,
                     'log2_fold_change':[1.]*10+[0.]*20})
    out=overrepresentation(de,{'test':set(de.gene_symbol[:12])},min_size=2)
    assert out.iloc[0].universe_size==30
    assert out.iloc[0].overlap_size==10
    workflow(tmp_path)
    assert (tmp_path/'01_workflow.png').stat().st_size>1000

def test_reactome_v97_and_legacy_formats(tmp_path):
    fixture=tmp_path/'sets.gmt'
    fixture.write_text('New\tR-HSA-1\tESR1\tPGR\nOld\tR-HSA-2\tReactome Pathway\tMKI67\tCCND1\n')
    sets=read_reactome(fixture)
    assert sets['New [R-HSA-1]']=={'ESR1','PGR'}
    assert sets['Old [R-HSA-2]']=={'MKI67','CCND1'}

def test_marker_absence_is_explicit():
    markers=pd.read_csv('results/marker_summary.csv').set_index('gene_symbol')
    assert len(markers)==18
    assert markers.loc['CYP19A1','measurement_status']=='not in processed matrix'
    assert pd.isna(markers.loc['CYP19A1','adjusted_p_value'])
