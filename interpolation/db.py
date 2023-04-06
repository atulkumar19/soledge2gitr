import sqlite3
import argparse
import os
import csv
import numpy as np
import collections
import copy
import warnings


def getSQLArrGenString(fName, dType, length):
    tString = ""
    if dType == int:
        tString = "INT"
    elif dType == float:
        tString = "REAL"
    else:
        raise Exception('Requested Incompatible SQL Trype')
    retStr = ""
    return retStr
    
    
def ProcessSOLEDGE(dbname, filepath,filename):

    fgDB = sqlite3.connect(dbname, timeout=45.0)
    fgCursor = fgDB.cursor()

    gndString = ""
    gndString = "CREATE TABLE IF NOT EXISTS SOLEDGE(R REAL, Z REAL, bphi REAL, br REAL, bz REAL, ve REAL, ne REAL, te REAL, vi REAL, ni REAL, ti REAL);" 
    

    sqlDB = sqlite3.connect(dbname)
    sqlCursor = sqlDB.cursor()
    sqlCursor.execute(gndString)

    txt = os.path.join(filepath, filename)

    trainingEntries = np.loadtxt(txt, comments='#')
    
    
    
    ZBARvalues = collections.namedtuple('Tablequantities', 'R Z bphi br bz ve ne te vi ni ti')
    for row in trainingEntries:
            res = ZBARvalues(R=row[0], Z=row[1], bphi=row[2],br=row[3],bz=row[4], ve=row[5], ne=row[6], te=row[7], vi=row[8], ni=row[9], ti=row[10])
            sqlCursor = sqlDB.cursor()
            insString = "INSERT INTO SOLEDGE VALUES(?, ?, ?,?, ?, ?, ?,?, ?, ?,?)"
            insArgs = (res.R, res.Z,res.bphi, res.br, res.bz, res.ve, res.ne,res.te,res.vi,res.ni,res.ti)
            sqlCursor.execute(insString, insArgs)
            sqlDB.commit()
            sqlCursor.close()

dbname='SOLEDGE.db'

B_path='.'   #current folder
B_filename='data.txt'  #soledge data


ProcessSOLEDGE(dbname, B_path,B_filename)
