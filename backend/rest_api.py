#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import sqlite3
from flask import Flask, request, redirect, url_for, send_from_directory, jsonify,g
from flask_restful import Resource, Api
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit import Chem
from rdkit.Chem import AllChem
from pythonds.basic.stack import Stack
import os
import json
import datetime
import time
from werkzeug.utils import secure_filename
from flask_cors import CORS
from database import *
UPLOAD_FOLDER = '/home/tanapat_ruengsatra/htdocs/input'
OUTPUT_FOLDER = '/home/tanapat_ruengsatra/htdocs/output'
ALLOWED_EXTENSIONS = set(['smi', 'pbd', 'sdf'])
OPTIONS = {'3d':0, '4d':1}


app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['OUTPUT_FOLDER'] = OUTPUT_FOLDER
app.config['CORS_HEADERS'] = 'Content-Type'
cors = CORS(app, resources={r"/*": {"origins": "*"}}, headers='Content-Type')

class ReactionItem:
    def __init__(self, id, name, rxn_smart, description):
        try:
            self.rxn =  AllChem.ReactionFromSmarts(str(rxn_smart))
        except Exception:
            self.rxn = None
        self.smart = rxn_smart
        self.id = id
        self.name = name
        self.description = description


def getAllReaction():
    rxn_list = []
    with app.app_context():
        all_reaction = get_all_reaction()
        for rxn in all_reaction:
            rxn_list.append(ReactionItem(rxn[0],rxn[1], rxn[2],rxn[3]))
    return rxn_list


def findReactionById(id):
    all_rxn_list = getAllReaction()
    for rxn in all_rxn_list:
        if(rxn.id == int(id)):
            print("match" + str(rxn.id))
            return rxn
    return None

class StackItem:
    def __init__(self, mol,num_reaction):
        self.mol = mol
        self.num_reaction = num_reaction

def runReaction(rxn, mol_reactant, reagent=None):
    products = []
    try:
        if reagent!=None:
            result = rxn.RunReactants((mol_reactant,reagent))
        else:
            result = rxn.RunReactants((mol_reactant,))
        try: 
            if result[0]:
                for product in result[0]:
                    products.append(product)       
        except IndexError:
            products.append(mol_reactant)
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
def make3D(mol_list, filename):
    remover = SaltRemover()
    outf = Chem.SDWriter(app.config['OUTPUT_FOLDER']+'/' + filename)
    for mol_index, mol in enumerate(mol_list):
        try:
            try:
                non_salt = remover(mol)
            except :
                non_salt = mol
            mol_hs = Chem.AddHs(non_salt)
            AllChem.EmbedMolecule(mol_hs, useRandomCoords=True)
            AllChem.MMFFOptimizeMolecule(mol_hs)
            outf.write(mol_hs)
        except:
            print('ERROR :'+ str(mol_index))
def readFile(filename):
    print('reading file')
    suppl = Chem.SmilesMolSupplier(filename,delimiter='\t',titleLine=False)
    mol_list = []
    for mol in suppl:
        mol_list.append(mol)
    print('Num of mol: '+ str(len(mol_list)))
    return mol_list
       
def executeFile(filename, input_rxn_list, option_list = []):
    mol_list = readFile(filename)
    for rxn in input_rxn_list:
        rxn_id = int(rxn['id'])
        if rxn.has_key('reagent'):
            print(str(rxn['reagent']))
            reagent = Chem.MolFromSmiles(str(rxn['reagent']))
        else:
            reagent = None
        rxnItem = findReactionById(rxn_id)
        rxn_object = rxnItem.rxn
        loop = int(rxn['loop'])
        print('Loop: '+str(loop))
        mol_list = runReactionList(rxn_object, mol_list, loop_num=loop, reagent=reagent)
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d-%H-%M-%S')
    output_file = 'output'+st+'.smi'
    created_output = False
    sample = ''
    if len(option_list)>0:
        for option in option_list:
            if option['id'] == OPTIONS['3d']:
                sdf_file = 'output'+st+'.sdf'
                print('3D output: '+sdf_file)
                result = sdf_file
                created_output = True
                for index, mol in enumerate(mol_list):
                    smile = Chem.MolToSmiles(mol)
                    if (index < 50):
                        sample += "%s compound%d\n"  %(smile, (index+1))
                make3D(mol_list, sdf_file)
               
    if(not created_output):
        w = open(app.config['OUTPUT_FOLDER']+'/'+output_file,"w+")
        #w = Chem.SmilesWriter(app.config['OUTPUT_FOLDER']+'/'+output_file)
        result = output_file
        for index, mol in enumerate(mol_list):
            smile = Chem.MolToSmiles(mol)
            w.write("%s compound%d\n"  %(smile, (index+1)))
            if (index < 50):
                sample += "%s compound%d\n"  %(smile, (index+1))
            #print(str(smile) + " compound " + str(index+1))
    return result, sample

@app.teardown_appcontext
def close_connection(exception):
    db = getattr(g, '_database', None)
    if db is not None:
        db.close()

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
        if 'file' not in request.files:
            return str({'status':'error'})  
        file = request.files['file']
        reaction = json.loads(request.form['reaction'])
        options = json.loads(request.form['options'])
        print(options)
        if file.filename == '':
            return json.dumps({'status':'no file'})
          
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(file_path)
            result, sample = executeFile(file_path, reaction, option_list=options)
            return json.dumps({"output":result, "sample": sample})


@app.route('/download/<path:filename>')
def download_file(filename):
    return send_from_directory(app.config['OUTPUT_FOLDER'],
                               filename, as_attachment=True)

@app.route('/resource/reaction', defaults={'rid': None}, methods=['GET','POST'] )
@app.route('/resource/reaction/<rid>', methods=['GET','POST', 'DELETE'])
def get_reaction(rid):
    if request.method == 'GET':
        all_rxn_list = getAllReaction()
        reaction_list = []
        for rxn in all_rxn_list:
            if(rid == None ):
                reaction_list.append({"id":rxn.id,"name":rxn.name, "smart":rxn.smart, "description": rxn.description})
            elif(int(rxn.id)==int(rid)):
                reaction_list.append({"id":rxn.id,"name":rxn.name, "smart":rxn.smart, "description": rxn.description})
        return jsonify(reaction=reaction_list)
    elif request.method == 'POST':
        reaction = json.loads(request.form['reaction'])
        name = reaction['name']
        smart = reaction['smart']
        description = reaction['description']
        with app.app_context():
            if(rid == None ):
                insert_db(name,smart,description)
            else:
                update_db(rid, name, smart, description)
        return json.dumps({"result":"insert successfully"})
    elif request.method == 'DELETE':
        print('delete ' + str(rid))
        with app.app_context():
            delete_db(rid)
        return json.dumps({"result":"delete successfully"})
    return json.dumps({"result":"please check http method"})
@app.route('/tools/smart/<smart>')
def check_smart(smart):
    try:
        rxn =  AllChem.ReactionFromSmarts(str(smart))
        return json.dumps({"result":True})
    except Exception:
        return json.dumps({"result":False})
if __name__ == "__main__":
    app.run(debug=True, host='0.0.0.0') 