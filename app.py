import os
from os.path import join
import sys
print(sys.version)
import numpy as np
import flask

from tqdm import tqdm
import dash
print(dash.__file__)
import dash_core_components as dcc
import pandas as pd
import dash_html_components as html
from dash.dependencies import Input, Output
from base64 import b64encode
from rdkit.Chem import Draw
from rdkit import Chem

def PrintAsBase64PNGString(mol, highlightAtoms=[], molSize = (200, 200)):
    data = Draw._moltoimg(mol,molSize,highlightAtoms,"",returnPNG=True, kekulize=True)
    return b64encode(data).decode('ascii')

def nan_NA(v):
    if v=='N/A':
        return np.nan
    else:
        return float(v)
    
app = dash.Dash(__name__)
server = app.server

mols = [m for m in Chem.SDMolSupplier('chembl_full.sdf', removeHs=False)]
df = pd.DataFrame([m.GetPropsAsDict() for m in tqdm(mols)])
df=df.assign(mol=mols)
df.CB1_pKi = df.CB1_pKi.apply(nan_NA)
df.CB2_pKi = df.CB2_pKi.apply(nan_NA)
df=df.dropna(subset=('CB1_pKi', 'CB2_pKi'))
print('Rendering structures...')
df=df.assign(imgs=[PrintAsBase64PNGString(m) for m in tqdm(df.mol.values)])

sliders = {}
for r in ['CB1', 'CB2']:
    sliders[r] = dcc.RangeSlider(id='%s-slider'%r,
                             min=df['%s_pKi'%r].min(),
                             max=df['%s_pKi'%r].max(),
                             step=0.1,
                             value=[df['%s_pKi'%r].min(),df['%s_pKi'%r].max()],
                             marks={n:{'label':'%.2f'%n} for n in np.linspace(df['%s_pKi'%r].min(), 
                                    df['%s_pKi'%r].max(), 5)})

app.layout = html.Div([
    html.Table([
            html.Tr([html.Td('CB1 pKi range:'), html.Td('CB2 pKi range:')]),
            html.Tr([
                    html.Td(sliders['CB1'], style={'width':'25%','padding-right':'10em'}),
                    html.Td(sliders['CB2'], style={'width':'25%','padding-right':'10em'})
                    ]),
            html.Tr([html.Td(html.Button('Show', id='button', style={'width':'15%','margin-top':'5em'})),
                     html.Td()])
    ] ,style={'width':'100%'}),
    dcc.Loading(html.Table([], 
                id='results',
                style={'width':'100%',
                       'border-spacing':'1em', 
                       'margin-top':'2em',
                       'border-top-width':'1px',
                       'border-top-style':'dashed',
                       'border-top-color':'grey'}))
    ])
   

    
@app.callback(
dash.dependencies.Output('results', 'children'),
[dash.dependencies.Input('button', 'n_clicks')],
[dash.dependencies.State('CB1-slider', 'value'),
 dash.dependencies.State('CB2-slider', 'value')])

def update_table(n_clicks, cb1, cb2):
    if not n_clicks:
        return []
    cb1_min, cb1_max = cb1
    cb2_min, cb2_max = cb2
    mask = (df.CB1_pKi>=cb1_min) & (df.CB1_pKi<=cb1_max)
    mask &= (df.CB2_pKi>=cb2_min) & (df.CB2_pKi<=cb2_max)
    imgs = [html.Img(src='data:image/png;base64,%s'%img) for img in df[mask].imgs]
    labels = ['CB1 pKi: %.2f, CB2 pKi: %.2f'%(cb1,cb2) for cb1,cb2 in zip(df.CB1_pKi.values,
              df.CB2_pKi.values)]
    tds = [html.Td([html.Tr(i), html.Tr(l)]) for i,l in zip(imgs, labels)]
    inds = np.array_split(np.arange(len(tds)), int(np.ceil(len(tds)/6)))
    trs = [html.Tr([tds[i] for i in ind]) for ind in inds]    
    
    return trs

    

if __name__ == '__main__':
    app.run_server(debug=False)