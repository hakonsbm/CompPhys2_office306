#Turbomole converter 
#Converts basis files in turbomole format to C++ files readable by Fermion Mingle

import glob
import os

def extract_m_basisinfo(filename, printing = False):
    #Extract multiple basises from a turbomole file
    swtch = 0
    f = open(filename)
    
    comments = []
    basis_sets = [] #list containing all basis sets
    basis = []
    basisW = []
    basisE = []
    basisT = []
    contractedtype = []
    activecontracted = -1
    for i in f:
        I = i.split()
        if len(I)==0:
            continue
        if I[0] == "*":
            continue #ignore wildcards
        if I[0] == '#':        
            comments.append("//")
            comments[-1] += i[1:]
            #for e in I:
            #    comments[-1] += e + " "
            continue
        #print I
        
        if I[0][0] == "$":
            continue
        if I[0] in ["h", "he", "li", "be", "b", "c", "n", "o", "f", "ne", "na" ,"mg", "al", "si", "p", "s", "cl", "ar", "k", "ca", "sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu", "zn", "ga", "ge", "as", "se","br", "kr", "rb", "sr", "y", "zr", "nb", "mo", "tc", "ru", "rh", "pd", "ag", "cd", "in", "sn", "sb", "te", "i", "xe"]:
            #If not first time here, put current basis into list
            if swtch == 1:
                basis_sets.append([bT, basisW, basisE, contractedtype, len(contractedtype)])
                basisW = []
                basisE = []
                contractedtype = []
                activecontracted = -1
            swtch = 1

            #Begin collecting data
            bT = I[1]+"_"+I[0]
            bT = bT.replace("-", "_") #basistype
            bT = bT.replace("(", "") #basistype
            bT = bT.replace(")", "") #basistype                        
            bT = bT.replace(",", "") #basistype                                    
            continue
        if I[0] in "123456789":
            if I[1] == "s":
                #print "Identified an s orbital."
                contractedtype.append(0)
            if I[1] == "p":
                #print "Identified an p orbital."
                contractedtype.append(1)
            if I[1] == "d":
                #print "Identified an d orbital."
                contractedtype.append(2)
            if I[1] == "f":
                #print "Identified an f orbital."
                contractedtype.append(3) 
            #Create contracted
            basis.append([])
            basisE.append([])
            basisW.append([])
            activecontracted += 1   
            #print I
        else:
            try:
                e = float(I[0].replace("D", "E"))
                w = float(I[1].replace("D", "E"))
                basisE[activecontracted].append(e)
                basisW[activecontracted].append(w)
                #basisE[activecontracted].append(float(I[1]))
                #basisW[activecontracted].append(float(I[0]))
            except:
                comments.append("//Ignored the following information from file:")
                for e in I:
                    comments[-1] += e + " "
                #comments[-1] += "Ignored the following information:" + I
                #print "Failed to read weigths and exponents."
                #print I
        #print I   
    if printing:
        print comments
        print contractedtype
        print len(contractedtype)
        print basisE
        print basisW
    f.close()
    basis_sets.append([bT, basisW, basisE, contractedtype, len(contractedtype)])
    #return comments, contractedtype, len(contractedtype), basisE, basisW
    return comments, basis_sets

def extract_basisinfo(filename, printing = False):
    swtch = 0
    f = open(filename)
    
    comments = []
    basis = []
    basisW = []
    basisE = []
    basisT = []
    contractedtype = []
    activecontracted = -1
    for i in f:
        I = i.split()
        if len(I)==0:
            continue
        if I[0] == "*":
            continue #ignore wildcards
        if I[0] == '#':        
            comments.append("//")
            for e in I:
                comments[-1] += e + " "
            continue
        #print I
        
        if I[0][0] == "$":
            swtch += 1
            swtch = swtch %2
            #if swtch == 1:
            #Add new function   
            continue
        if I[0] in ["h", "he", "li", "be", "b", "c", "n", "o", "f", "ne", "na" ,"mg", "al", "si", "p", "s", "cl", "ar", "k", "ca", "sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu", "zn", "ga", "ge", "as", "se","br", "kr", "rb", "sr", "y", "zr", "nb", "mo", "tc", "ru", "rh", "pd", "ag", "cd", "in", "sn", "sb", "te", "i", "xe"]:
            bT = I[1]+"_"+I[0]
            bT = bT.replace("-", "_")
            basisT.append(bT)
            print basisT[-1]
            continue
        if I[0] in "123456789":
            if I[1] == "s":
                #print "Identified an s orbital."
                contractedtype.append(0)
            if I[1] == "p":
                #print "Identified an p orbital."
                contractedtype.append(1)
            if I[1] == "d":
                #print "Identified an d orbital."
                contractedtype.append(2)
            if I[1] == "f":
                #print "Identified an f orbital."
                contractedtype.append(3) 
            #Create contracted
            basis.append([])
            basisE.append([])
            basisW.append([])
            activecontracted += 1   
            #print I
        else:
            try:
                w = float(I[0].replace("D", "E"))
                e = float(I[1].replace("D", "E"))
                basisE[activecontracted].append(e)
                basisW[activecontracted].append(w)
            except:
                comments.append("//Ignored the following information from file:")
                for e in I:
                    comments[-1] += e + " "
                #comments[-1] += "Ignored the following information:" + I
                #print "Failed to read weigths and exponents."
                #print I
        #print I          
    if printing:
        print comments
        print contractedtype
        print len(contractedtype)
        print basisE
        print basisW
    f.close()
    return comments, contractedtype, len(contractedtype), basisE, basisW



def createfunction(fname, comments, orbital_types, N_orbitals, exponents, weights, coeffs):
    endline = ["\n"]
    cppclass  = []
    cppheader = []
    for i in comments:
        cppclass.append(i[:-2])
        #cppheader.append(i)        
        
    #cppclass.append(endline)
    #cppheader.append(endline)
    
    #writing to header
    cppheader.append("    void add_"+fname+"(const vec corePos, const vec c);")

    

    cppclass.append("void basisbank::add_"+fname+"(const vec corePos, const vec c){")
    orb_count = 0
    for i in range(N_orbitals):
        if orbital_types[i] == 0:
            #create an s-orbital
            cppclass.append("    // s-orbital")
            #cppclass.append("    bs.add_state();")
            
            for e in range(len(exponents[i])):
                cppclass.append("    contracted()->add_primitive(%.8f,%.8f,0,0,0,corePos);" % (exponents[i][e], weights[i][e]))
                orb_count += 1
            cppclass.append("    contracted()->contract_orb(c[%d]);" %i)        
            
            """
            if len(exponents[i]) == 3:
                cppclass.append("    contracted()->contract_orb_2s();")
                orb_count = 0
            elif len(exponents[i]) == 2:
                cppclass.append("    contracted()->contract_orbK();")
            else:
                cppclass.append("    contracted()->contract_orb_1s();")
            """
                

            

        if orbital_types[i] == 1:
            #create an p-orbital
            for coor in [["1,0,0", 0],["0,1,0", 1],["0,0,1", 2]]:
                cppclass.append("    // p-orbital")
                for e in range(len(exponents[i])):
                    cppclass.append("    contracted()->add_primitive(%.8f,%.8f,%s,corePos);" % (exponents[i][e], weights[i][e],coor[0]))
                cppclass.append("    contracted()->contract_orb(c[%d]);" %(i+2*int(coor[1]))) 

            
        """
                if len(exponents[i]) == 1:
                    cppclass.append("    contracted()->contract_orb();")
                else:
                    cppclass.append("    contracted()->contract_orbK();")
            """
        

        ######## NOT IMPLEMENTED ########
        """
        if orbital_types[i] == 2:
            #create an d-orbital
            cppclass.append("    // d-orbital (2)")
            cppclass.append("    bs.add_state();")
            for e in range(len(exponents[i])):
                cppclass.append("    contracted()->add_primitive(%.8f,%.8f,2,0,0,corePos);" % (exponents[i][e], weights[i][e]))
                cppclass.append("    contracted()->add_primitive(%.8f,%.8f,0,2,0,corePos);" % (exponents[i][e], weights[i][e]))
                cppclass.append("    contracted()->add_primitive(%.8f,%.8f,0,0,2,corePos);" % (exponents[i][e], weights[i][e]))

                cppclass.append("    contracted()->add_primitive(%.8f,%.8f,1,1,0,corePos);" % (exponents[i][e], weights[i][e]))
                cppclass.append("    contracted()->add_primitive(%.8f,%.8f,0,1,1,corePos);" % (exponents[i][e], weights[i][e]))
                cppclass.append("    contracted()->add_primitive(%.8f,%.8f,1,0,1,corePos);" % (exponents[i][e], weights[i][e]))

        if orbital_types[i] == 3:
            #create an d-orbital
            cppclass.append("    // d-orbital (3)")
            cppclass.append("    bs.add_state();")
            for e in range(len(exponents[i])):
                cppclass.append("    contracted()->add_primitive(%.8f,%.8f,3,0,0,corePos);" % (exponents[i][e], weights[i][e]))
                cppclass.append("    contracted()->add_primitive(%.8f,%.8f,0,3,0,corePos);" % (exponents[i][e], weights[i][e]))
                cppclass.append("    contracted()->add_primitive(%.8f,%.8f,0,0,3,corePos);" % (exponents[i][e], weights[i][e]))

                cppclass.append("    contracted()->add_primitive(%.8f,%.8f,1,2,0,corePos);" % (exponents[i][e], weights[i][e]))
                cppclass.append("    contracted()->add_primitive(%.8f,%.8f,0,1,2,corePos);" % (exponents[i][e], weights[i][e]))
                cppclass.append("    contracted()->add_primitive(%.8f,%.8f,1,0,2,corePos);" % (exponents[i][e], weights[i][e]))

                cppclass.append("    contracted()->add_primitive(%.8f,%.8f,2,1,0,corePos);" % (exponents[i][e], weights[i][e]))
                cppclass.append("    contracted()->add_primitive(%.8f,%.8f,0,2,1,corePos);" % (exponents[i][e], weights[i][e]))
                cppclass.append("    contracted()->add_primitive(%.8f,%.8f,2,0,1,corePos);" % (exponents[i][e], weights[i][e]))

                cppclass.append("    contracted()->add_primitive(%.8f,%.8f,1,1,1,corePos);" % (exponents[i][e], weights[i][e]))
        """

    cppclass.append("}\n")
    
    return cppheader, cppclass

def saveclass(cppclasses, cppheaders):
    cppfile = "//This file is maintained by an external python script and should not be edited manually.\n\
#include <armadillo>\n\
#include \"basisbank.h\"\n\
#include \"contracted.h\"\n\n\
using namespace std;\n\
using namespace arma;\n\n\
basisbank::basisbank()\n{\n    \
initContracted(new Contracted);\n\
}\n\
\n\
\n"
    
    for i in cppclasses:
        for e in i:
            cppfile += e+ "\n"
        cppfile += "\n"
    cppfile += "double basisbank::get_orb() {return contracted()->get_orb();}\n"
    """
    cppfile += "double basisbank::get_orb_1s() {return contracted()->get_orb_1s();}\n"
    cppfile += "double basisbank::get_orb_2s() {return contracted()->get_orb_2s();}\n"
    cppfile += "double basisbank::get_orb_pX() {return contracted()->get_orb_pX();}\n"
    cppfile += "double basisbank::get_orb_pY() {return contracted()->get_orb_pY();}\n"
    cppfile += "double basisbank::get_orb_pZ() {return contracted()->get_orb_pZ();}\n"
    """
    #print cppfile
    
    
    cppheader = "//This file is maintained by an external python script and should not be edited manually.\n\
#ifndef BASISBANK_H\n\
#define BASISBANK_H\n\
#include <armadillo>\n\
#include \"contracted.h\"\n\n\
using namespace std;\n\
using namespace arma;\n \n\
class basisbank{\n\
public:\n"
    cppheader += "    basisbank();\n"
    for i in cppheaders:
        for e in i:
            cppheader += e + "\n"
    cppheader += "    void initContracted(Contracted *contracted) {m_contracted = contracted;}\n"
    cppheader += "    Contracted *contracted() {return m_contracted;}\n"
    cppheader += "    double get_orb();\n"
    """
    cppheader += "    double get_orb_1s();\n"
    cppheader += "    double get_orb_2s();\n"
    cppheader += "    double get_orb_pX();\n"
    cppheader += "    double get_orb_pY();\n"
    cppheader += "    double get_orb_pZ();\n"

    """

    cppheader += "private:\n"
    cppheader += "    Contracted *m_contracted;\n"

    cppheader += "};\n"
    cppheader += "#endif // BASISBANK_H"
    #print cppheader
    #Write files
    cpp_header = open("basisbank.h", "w")
    cpp_class = open("basisbank.cpp", "w")
    cpp_header.write(cppheader)
    cpp_class.write(cppfile)
    cpp_header.close()
    cpp_class.close()

    
def convertfolder():
    
    cppclasses = []
    cppheaders = []
    folder  = [] #a list containing all files in folder
    #os.chdir("/")
    #print os.getcwd()
    for filename in os.listdir(os.getcwd()):
        if filename.endswith(".txt"):
            #try:
            #comments, orbital_types, N_orbitals, exponents, weigths = extract_basisinfo(filename)
            comments, I = extract_m_basisinfo(filename, False)
            commented = 1
            #print comments
            #print "I:", I

            for i in I:
                #i = [bT, basisW, basisE, contractedtype, len(contractedtype)]
                if commented == 0:
                    cppheader, cppclass = createfunction(i[0], comments, i[3], i[4], i[2], i[1], [])
                    cppclasses.append(cppclass)
                    cppheaders.append(cppheader)
                    commented = 1
                else:
                    cppheader, cppclass = createfunction(i[0], [""], i[3], i[4], i[2], i[1], [])
                    cppclasses.append(cppclass)
                    cppheaders.append(cppheader)

                    
            #except:
                    #print "Failed to import file:", filename

    saveclass(cppclasses, cppheaders)


convertfolder() #converts all .txt files in current folder