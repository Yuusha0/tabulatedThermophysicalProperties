#!/usr/bin/python3
# -*- coding: utf-8 -*-

from re import sub

class thermophysicalTable() :
    """Base class for manipulate thermophysical table
    needed by OpenFOAM tabulated thermophysical properties model."""

    def __init__(self) :
        self.table = []
        self.fileName = ''

    def read(self, fileName='') :
        """ Read a tabulated thermophysical file and store it into table"""

        with open(fileName, 'r') as f :
            for line in f :
                tmpList = sub('\(+|\)+', ' ', line).strip().split()
                try :
                    count = len(tmpList) - 1
                    tupleList=[]
                    for i in range(1,len(tmpList)-1,2):
                        tupleList.append((tmpList[i], tmpList[i+1]))
                    self.table.append([tmpList[0], tupleList])
                except IndexError :
                    continue

    def write(self, fileName='') :
        """ Write a thermopysical file """

        with open(fileName, 'w') as f :
            f.write('(\n')
            for element in self.table :
                stringOut=""
                for elem in element[1][:] :
                    stringOut += '({} {}) '.format(elem[0], elem[1])
                f.write('( {} ({}))\n'.format(element[0], stringOut))
            f.write(')')

    def importation(self, fileName='', ext = '') :
        """ Import thermophysical data from a file.
        Only CSV is implemented yet """

    def transpose(self) :
        """ Invert lines and columns in a tabulated thermophysical list """

        tableTmp = list(self.table)
        self.table = []
        for element in tableTmp :
            for elem in element[1] :
                find = False
                for pres in self.table :
                    if elem[0] == pres[0] :
                        pres[1].append((element[0], elem[1]))
                        find = True
                        break
                if not find :
                    self.table.append([elem[0], [(element[0], elem[1])]])


thermo = thermophysicalTable()
thermo.read('hTable_orig')
thermo.transpose()
thermo.write('hTable_transposed')
