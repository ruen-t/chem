import sqlite3
from flask import g,Flask

DATABASE = '/home/tanapat_ruengsatra/htdocs/Chem.db'
def get_db():
    db = getattr(g, '_database', None)
    if db is None:
        db = g._database = sqlite3.connect(DATABASE)
    return db
def query_db(query, args=(), one=False):
    cur = get_db().execute(query, args)
    rv = cur.fetchall()
    cur.close()
    return (rv[0] if rv else None) if one else rv
def insert_db(name, smart, description):
    db = get_db()
    cur = db.cursor()
    cur.execute("INSERT INTO Reaction  VALUES (NULL,?,?,?)",(name,smart,description))
    db.commit()
    cur.close()
def update_db(key, name, smart, description):
    db = get_db()
    cur = db.cursor()
    cur.execute("UPDATE Reaction SET name=?, smart=?, description=? WHERE id=?",(name,smart,description,key))
    db.commit()
    cur.close()
def delete_db(key):
    db = get_db()
    cur = db.cursor()
    cur.execute("DELETE FROM Reaction WHERE id=?",(str(key)))
    db.commit()
    cur.close()
def get_all_reaction():
     return query_db("select * from Reaction")