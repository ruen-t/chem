#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import sqlite3
from flask import Flask, request, redirect, url_for, send_from_directory, jsonify,g
from flask_restful import Resource, Api
from rdkit import Chem
import os
import json
import datetime
import time
from werkzeug.utils import secure_filename
from flask_cors import CORS
from database import *
from ChemProcess import *
UPLOAD_FOLDER = '/home/tanapat_ruengsatra/htdocs/input'
OUTPUT_FOLDER = '/home/tanapat_ruengsatra/htdocs/output'
ALLOWED_EXTENSIONS = set(['smi', 'pbd', 'sdf'])
OPTIONS = {'3d':0, 'prepare':1}
SAMPLE_PREFIX = 'sample_'


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

def cleanMolList(mol_list):
    cleaned_mol_list = []
    for index, mol in enumerate(mol_list):
        try:
            name = mol.GetProp("_Name")
        except:
            name = "unspecified"+str(index)
            mol.SetProp("_Name", name)
        cleaned_mol_list.append(mol)
    return cleaned_mol_list

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

def readFile(filename):
    smiles_file = open(filename,"r")
    mol_list = []
    for index, line in enumerate(smiles_file):
        split_text = line.split('\t')
        if(len(split_text) < 2):
            split_text = line.split(' ')
            if(len(split_text) < 2):
                split_text = [line]
        smiles = split_text[0]
        name = ''
        for i in range(len(split_text)):
            if(i>0):
                name += split_text[i]
        name = name.rstrip()
        if len(name) == 0:
            name = "unspecified_"+str(index)
        mol = Chem.MolFromSmiles(smiles)
        mol.SetProp('_Name', name)
        mol_list.append(mol)
    print('Num of mol: '+ str(len(mol_list)))
    return mol_list


def make3DAfterReaction(mol_list, sdf_file, option, sample_size=50):
    filename = app.config['OUTPUT_FOLDER']+'/'+sdf_file
    smile_list, name_list = makeSmileList(mol_list)
    new_mol_list = []
    sample_mol_list = []
    sample_filename =  app.config['OUTPUT_FOLDER']+'/'+SAMPLE_PREFIX+sdf_file
    for i  in range(len(smile_list)):
        new_mol = Chem.MolFromSmiles(smile_list[i])
        new_mol.SetProp("_Name", name_list[i])
        new_mol.SetProp("Name", name_list[i])
        new_mol_list.append(new_mol)
        if(i < sample_size):
            sample_mol_list.append(new_mol)
    removeSalt = option['removeSalt']
    ionize = option['ionize']
    pH = option['pH']
    make3D(sample_mol_list, sample_filename, removeSalt, ionize=ionize, pH=pH)
    return make3D(new_mol_list, filename, removeSalt, ionize=ionize, pH=pH), (SAMPLE_PREFIX+sdf_file)


def writeMolToSmileFile(mol_list, filename):
    sample = ''
    w = open(app.config['OUTPUT_FOLDER']+'/'+filename,"w+")
    for index, mol in enumerate(mol_list):
        smile = Chem.MolToSmiles(mol)
        w.write("%s\t%s\n"  %(smile, mol.GetProp("_Name")))
        if (index < 50):
            sample += "%s\t%s\n"  %(smile, mol.GetProp("_Name"))
    return filename, sample

def getSample(mol_list, number):
    sample = ''
    for index, mol in enumerate(mol_list):
        smile = Chem.MolToSmiles(mol)
        if (index < number):
            sample += "%s\t%s\n"  %(smile, mol.GetProp("_Name"))
        else:
            break
    return sample
def runOption(mol_list, option):
    outputCreated = False
    ts = time.time()
    sdf_file = None
    sample = ""
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d-%H-%M-%S')
    if option['id'] == OPTIONS['3d']:
        sdf_file ='output'+st+'.sdf'
        mol_list_result, sample = make3DAfterReaction(mol_list, sdf_file, option)
        outputCreated = True
        #sample = getSample(mol_list_result, 50)
    if option['id'] == OPTIONS['prepare']:
        ionize = option['ionize']
        pH = option['pH']
        print("pH: "+str(pH))
        addHs = option['addHs']
        mol_list_result = prepareMol(mol_list, ionize=ionize, pH = float(pH), addHs = addHs)

    return sdf_file, sample, outputCreated, mol_list_result
def executeFile(filename, input_rxn_list, option_list = [], cleanFile = True):
    mol_list = readFile(filename)
    sdf_output = False
    if cleanFile:
        mol_list = cleanMolList(mol_list)
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
    outputCreated = False
    if len(option_list)>0:
        for option in option_list:
            file_name, sample, outputCreated, mol_list = runOption(mol_list, option)
            if(outputCreated):
                sdf_output = True
    
    if (not sdf_output):
        file_name, sample = writeMolToSmileFile(mol_list, output_file)      
    return file_name, sample, sdf_output

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
        cleanInput = request.args.get("cleanInput")
        print(options)
        if file.filename == '':
            return json.dumps({'status':'no file'})
          
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(file_path)
            result, sample, sdf_output = executeFile(file_path, reaction, option_list=options, cleanFile = False)
            return json.dumps({"output":result, "sample": sample, "sdf": sdf_output})


@app.route('/download/<path:filename>')
def download_file(filename):
    return send_from_directory(app.config['OUTPUT_FOLDER'],
                               filename, as_attachment=True)

@app.route('/viewfile/<path:filename>')
def view_file(filename):
    return send_from_directory(app.config['OUTPUT_FOLDER'],
                               filename, as_attachment=False, mimetype='text/plain')


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

@app.route('/tools/smarts')
def check_smart():
    try:
        smarts = request.args.get("q")
        print(smarts)
        rxn =  AllChem.ReactionFromSmarts(str(smarts))
        return json.dumps({"result":True})
    except Exception:
        return json.dumps({"result":False})

@app.route('/tools/3d/<smiles>')
def getSDF(smiles):
    remover = SaltRemover()
    outf = Chem.SDWriter(OUTPUT_FOLDER+'/test.sdf')
    mol = Chem.MolFromSmiles(smiles)
    make3DFromMol(mol, outf, remover)
    return send_from_directory(app.config['OUTPUT_FOLDER'],
                               'test.sdf', as_attachment=True)
if __name__ == "__main__":
    app.run(debug=True, host='0.0.0.0') 