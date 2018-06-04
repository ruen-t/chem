from rdkit.Chem.SaltRemover import SaltRemover
from vs_utils.features import MolPreparator
from rdkit import Chem
from rdkit.Chem import AllChem
from pythonds.basic.stack import Stack
class StackItem:
    def __init__(self, mol,num_reaction):
        self.mol = mol
        self.num_reaction = num_reaction

def runReaction(rxn, mol_reactant, reagent=None):
    products = []
    name = mol_reactant.GetProp("_Name")
    try:
        if reagent!=None:
            result = rxn.RunReactants((mol_reactant,reagent))
        else:
            result = rxn.RunReactants((mol_reactant,))
        try: 
            if result[0]:
                for product in result[0]:
                    product.SetProp("_Name", name+"*")
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

def make3DFromMol(mol, sdWriter, saltRemover = None):
    if saltRemover!= None:
        try:
            non_salt = saltRemover(mol)
        except:
            non_salt = mol
    try:
        mol_hs = Chem.AddHs(non_salt)
    except:
        mol_hs = non_salt
    try:
        AllChem.EmbedMolecule(mol_hs, useRandomCoords=True)
    except:
        pass
    try:
        AllChem.MMFFOptimizeMolecule(mol_hs)
    except:
        pass
    sdWriter.write(mol_hs)
    return mol_hs

def make3D(mol_list, filename, removeSalt=True, ionize=False, pH = 7.0):
    mol_list_result = []
    mol_list = prepareMol(mol_list, ionize, pH, False )
    if removeSalt:
        remover = SaltRemover()
    else:
        remover = None
    outf = Chem.SDWriter(filename)
    for mol_index, mol in enumerate(mol_list):
        try:
            mol_result = make3DFromMol(mol, sdWriter = outf, saltRemover = remover)
            mol_list_result.append(mol_result)
            print("write mol")
        except:
            print('ERROR :'+ str(mol_index))
    return mol_list_result

def makeSmileList(mol_list):
    smile_list = []
    name_list = []
    for index, mol in enumerate(mol_list):
        name = mol.GetProp("_Name")
        smile = Chem.MolToSmiles(mol)
        smile_list.append(smile)
        name_list.append(name)
    return smile_list, name_list

def prepareMol(mol_list, ionize, pH, addHs ):
    prepare = MolPreparator(ionize=ionize, pH=pH, add_hydrogens=addHs)
    prepared_mol_list = []
    for mol in mol_list:
        try:
            prepared_mol = prepare(mol)
        except:
            prepared_mol = mol
        name = mol.GetProp("_Name")
        prepared_mol.SetProp("_Name", name)
        prepared_mol.SetProp("Name", name)
        prepared_mol_list.append(prepared_mol)
    return prepared_mol_list