#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from flask import Flask, request, redirect, url_for, send_from_directory, jsonify
from flask_restful import Resource, Api
from rdkit import Chem
from rdkit.Chem import AllChem
from pythonds.basic.stack import Stack
import os
import json
from werkzeug.utils import secure_filename
from flask_cors import CORS
UPLOAD_FOLDER = '/home/tanapat_ruengsatra/htdocs/input'
OUTPUT_FOLDER = '/home/tanapat_ruengsatra/htdocs/output'
ALLOWED_EXTENSIONS = set(['smi', 'pbd', 'jpg'])


app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['OUTPUT_FOLDER'] = OUTPUT_FOLDER
app.config['CORS_HEADERS'] = 'Content-Type'
cors = CORS(app, resources={r"/*": {"origins": "*"}}, headers='Content-Type')
# In[130]:
reaction = [
        {
            "id": 0,
            "name": "reaction0",
            "smart": "[O:2]=[C:1][O:3][C;H3]>>[O:2]=[C:1][O:3]",
            "reagent": None
        },
        {
            "id": 1,
            "name": "reaction1",
            "smart": "[O:2]=[C:1][O:3][C;H2][C;H3]>>[O:2]=[C:1][O:3]",
            "reagent": None
        },
        {
            "id": 2,
            "name": "reaction2",
            "smart": "[O:2]=[C:1][O:3][C;H2][C;H2][C;H3]>>[O:2]=[C:1][O:3]",
            "reagent": None
        },
        {
            "id": 3,
            "name": "reaction3",
            "smart": "[O:2]=[C:1][O:3][C;H2][C;H2][C;H2][C;H3]>>[O:2]=[C:1][O:3]",
            "reagent": None
        },
        {
            "id": 4,
            "name": "reaction4",
            "smart": "[O:2]=[C:1][OH].[N:3][C:4]>>[O:2]=[C:1][N:3][C:4]",
            "reagent": None
        }
    ]
class ReactionItem:
    def __init__(self, id, name, rxn_smart, reagent):
        self.rxn =  AllChem.ReactionFromSmarts(rxn_smart)
        self.id = id
        self.reagent = reagent
        self.name = name


def initReactionList():
    rxn_list = []
    for rxn in reaction:
        rxn_list.append(ReactionItem(rxn['id'],rxn['name'], rxn['smart'], rxn['reagent']))
    return rxn_list

rxn_list = initReactionList()
def findReactionById(id):
    for rxn in rxn_list:
        print("id:"+str(rxn.id))
        if(rxn.id == int(id)):
            print("match")
            return rxn
    return None

class StackItem:
    def __init__(self, mol,num_reaction):
        self.mol = mol
        self.num_reaction = num_reaction


# In[131]:

def runReaction(rxn, mol_reactant, reagent=None):
    products = []
    #print("run: " + Chem.MolToSmiles(mol_reactant))
    try:
        result = rxn.RunReactants((mol_reactant,reagent))
        try: 
            if result[0]:
                for product in result[0]:
                    #result_smile = Chem.MolToSmiles(product)
                    products.append(product)
            #print(str(products[0])+' compound'+ str(index+1))        
        except IndexError:
            products.append(mol_reactant)
            #print(str(smile)+' compound'+ str(index+1))
    except Exception:
        products.append(mol_reactant)
    return products



def loopReaction(rxn, mol_ractant,reagent = None, loop_num = 1):
    reactant_stack = Stack()
    all_products = []
    item = StackItem(mol = mol_ractant, num_reaction = 0)
    reactant_stack.push(item)
    while not reactant_stack.isEmpty():
        reactant = reactant_stack.pop()
        mol = reactant.mol
        num_reaction = reactant.num_reaction
        if(num_reaction == loop_num):
            all_products.append(mol)
        else:
            products = runReaction(rxn, mol,reagent)
            num_reaction = num_reaction+1
            for p in products:
                new_item = StackItem(mol = p, num_reaction = num_reaction)
                reactant_stack.push(new_item)
    return all_products



def runReactionList(rxn, mol_list,reagent = None, loop_num = 1):
    all_products = []
    for index, mol in enumerate(mol_list):
        products = loopReaction(rxn, mol, reagent, loop_num)
        for p in products:
            all_products.append(p)
    return all_products

def readFile(filename):
    print('reading file')
    suppl = Chem.SmilesMolSupplier(filename,delimiter='\t',titleLine=False)
    mol_list = []
    for mol in suppl:
        mol_list.append(mol)
    print('Num of mol: '+ str(len(mol_list)))
    return mol_list
        
def executeFile(filename, input_rxn_list):
    mol_list = readFile(filename)
    for rxn in input_rxn_list:
        rxn_id = rxn['reaction']
        if rxn.reagent:
            print(rxn['reagent'])
            reagent = Chem.MolFromSmiles(rxn['reagent'])
        else:
            reagent = None
        rxnItem = findReactionById(rxn_id)
        rxn = rxnItem.rxn
        mol_list = runReactionList(rxn, mol_list, loop_num=2, reagent=reagent)
    w = Chem.SmilesWriter(app.config['OUTPUT_FOLDER']+'/output.smi')
    result = app.config['OUTPUT_FOLDER']+'/output.smi'
    for index, mol in enumerate(mol_list):
        smile = Chem.MolToSmiles(mol)
        print(smile)
        w.write(mol)
        #result += str(smile)+' compound'+ str(index+1)
    return result


@app.route('/')
def hello_world():
    return str({'status':'upload successfully'})
    app.run(debug=True)
@app.route('/test')
def test():
    return str({'status':'test successfully'})  
    
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/upload', methods=['POST'])
def upload_file():
    if request.method == 'POST':
        # check if the post request has the file part
        if 'file' not in request.files:
            return str({'status':'error'})  
        file = request.files['file']
        reaction = request.form['reaction']
        for r in json.loads(reaction):
            print(r['id'])
        #reaction = request.args.get('reaction', None)
        #reaction_list = reaction.split(',')
        #for r in reaction_list:
        #    print(r)
        # if user does not select file, browser also
        # submit a empty part without filename
        if file.filename == '':
            return str({'status':'no file'})
          
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(file_path)
            result = executeFile(file_path, reaction)
            return str({'result':result})


@app.route('/download/<path:filename>')
def download_file(filename):
    return send_from_directory(app.config['OUTPUT_FOLDER'],
                               filename, as_attachment=True)

@app.route('/resource/reaction')
def get_reaction():
    return jsonify(reaction = reaction)
if __name__ == "__main__":
    app.run(debug=True, host='0.0.0.0') 