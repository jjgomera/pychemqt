#import sys
#from code import InteractiveConsole
#
#
#class FileCacher:
#    "Cache the stdout text so we can analyze it before returning it"
#    def __init__(self): self.reset()
#    def reset(self): self.out = []
#    def write(self,line): self.out.append(line)
#    def flush(self):
#        output = '\n'.join(self.out)
#        self.reset()
#        return output
#
#class Shell(InteractiveConsole):
#    "Wrapper around Python that can filter input/output to the shell"
#    def __init__(self):
#        self.stdout = sys.stdout
#        self.cache = FileCacher()
#        InteractiveConsole.__init__(self)
#        return
#
#    def get_output(self): sys.stdout = self.cache
#    def return_output(self): sys.stdout = self.stdout
#
#    def push(self,line):
#        self.get_output()
#        # you can filter input here by doing something like
#        # line = filter(line)
#        InteractiveConsole.push(self,line)
#        self.return_output()
#        output = self.cache.flush()
#        # you can filter the output here by doing something like
#        # output = filter(output)
#        print output # or do something else with it
#        return 
#
#if __name__ == '__main__':
#     sh = Shell()
#     sh.interact()

#import sqlite3
#from PyQt5 import QtGui
#
#connection = sqlite3.connect('elemental.db')
#cursor = connection.cursor()
#cursor.execute("""CREATE TABLE IF NOT EXISTS ELEMENTS \
#               (id INTEGER PRIMARY KEY,
#               name TEXT,
#               altname TEXT,
#               symbol TEXT,
#               serie TEXT,
#               group_element INTEGER,
#               period INTEGER,
#               block TEXT,
#               density_Solid TEXT,
#               density_Liq TEXT,
#               density_Gas TEXT,
#               appearance TEXT,
#               date TEXT,
#               country TEXT,
#               discover TEXT,
#               etymology TEXT,
#               atomic_mass TEXT,
#               atomic_volume TEXT,
#               atomic_radius TEXT,
#               covalent_radius TEXT,
#               vanderWaals_radius TEXT,
#               ionic_radii TEXT,
#               lattice_type TEXT,
#               space_group TEXT,
#               lattice_edges TEXT,
#               lattice_angles TEXT,
#               electron_configuration TEXT,
#               oxidation TEXT,
#               electronegativity TEXT,
#               electron_affinity TEXT,
#               first_ionization TEXT,
#               Tf TEXT,
#               Tb TEXT,
#               Heat_f TEXT,
#               Heat_b TEXT,
#               Cp TEXT,
#               k TEXT,
#               T_debye TEXT,
#               color TEXT,
#               notes TEXT)""")
#
#
#def removeSpaces(string):
#    while string[0] in ["\t", " "]:
#        string = string[1:]
#    while string[-1] in ["\t", " "]:
#        string = string[:-1]
#    return string
#
#def clean(string):
#    string = removeSpaces(string)
#    if string in ["Q_NA", "Q_UNK"]:
#        string = ""
#    if string[0:4] == 'N_("':
#        string = string[4:-2]
#    if string and string[0] == '"' and string[-1] == '"':
#        string = string[1:-1]
#    return string
#        
#def toNum(string, func=float):
#    if string:
#        try:
#            num = func(string)
#        except ValueError:
#            num = float(string.split("(")[1].split(",")[0])
#    else:
#        num = 0
#    return num
#    
#with open("elemnt.data", "r") as file:
#    for line in file:
#        element = {}
#        if line == "\n":
#            continue
#        data = line.split(",\t")
#        element["name"] = removeSpaces(data[0]).split('"')[1]
#        element["altname"] = clean(data[1])
#        element["symbol"] = data[2].split('"')[1]
#        id = int(removeSpaces(data[3]))
#        element["id"] = id
#        serie = clean(data[4])
#        if serie:
#            serie  = serie.split("::")[1]
#        element["serie"] = serie.capitalize().replace("_"," ")
#        group = toNum(clean(data[5]))
#        if not group:
#            group = 0
#        element["group_element"] = group
#        element["period"] = clean(data[6])
#        block = clean(data[7])
#        if block:
#            block  = block.split("::")[1]
#        element["block"] = block
#        
#        element["density_Solid"] = toNum(clean(data[8]))
#        element["density_Liq"] = toNum(clean(data[9]))
#        element["density_Gas"] = toNum(clean(data[10]))
#
#        element["appearance"] = clean(data[11])
#
#        event = clean(data[12])
#        if event:
#            date, country = event[7:-1].split(",")
#            country = clean(country)
#        else:
#            date, country = "", ""
#        element["date"] = date
#        element["country"] = country
#        element["discover"] = clean(data[13])
#        element["etymology"] = clean(data[14])
#        element["atomic_mass"] = toNum(clean(data[15]))
#        element["atomic_volume"] = toNum(clean(data[16]))
#        element["atomic_radius"] = toNum(clean(data[17]))
#        element["covalent_radius"] = toNum(clean(data[18]))
#        element["vanderWaals_radius"] = toNum(clean(data[19]))
#        element["ionic_radii"] = clean(data[20])
#        lattice = clean(data[21])
#        if lattice:
#            lattice  = lattice.split("::")[1]
#        element["lattice_type"] = lattice
#        element["space_group"] = toNum(clean(data[22]), int)
#        
#        txt = clean(data[23])
#        if txt:
#            edge = txt.split(",")[1:]
#            edge[-1] = edge[-1][:-1]
#            edges = [float(e) for e in edge]
#        element["lattice_edges"] = edges
#        
#        txt = clean(data[24])
#        if txt:
#            angle = txt.split(",")[1:]
#            angle[-1] = angle[-1][:-1]
#            angles = [float(e) for e in angle]
#        element["lattice_angles"] = angles
#        
#        configuration = clean(data[25])
#        if configuration and configuration[0] == '"':
#            configuration = configuration [1:-1]
#        element["electron_configuration"] = configuration
#        txt = clean(data[26])
#        if txt:
#            oxidations = txt.split(",")[1:]
#            oxidations[-1] = oxidations[-1][:-1]
#        element["oxidation"] = ", ".join(oxidations)
#        
#        element["electronegativity"] = toNum(clean(data[27]))
#        element["electron_affinity"] = toNum(clean(data[28]))
#        element["first_ionization"] = toNum(clean(data[29]))
#        element["Tf"] = toNum(clean(data[30]))
#        element["Tb"] = toNum(clean(data[31]))
#        element["Heat_f"] = toNum(clean(data[32]))
#        element["Heat_b"] = toNum(clean(data[33]))
#        element["Cp"] = toNum(clean(data[34]))
#        element["k"] = toNum(clean(data[35]))
#        element["T_debye"] = toNum(clean(data[36]))
#        color = clean(data[37])[7:-1].split(",")
#        try:
#            r, g, b = map(float, color)
#        except ValueError:
#            r, g, b = 1, 1, 1
#        color = QtGui.QColor(r*255, g*255, b*255)  
#        element["color"] = color.name()
#        element["notes"] = (clean(data[38]))
#
#        query = "insert into ELEMENTS ("
#        query2 = ") values ("
#        keys = []
#        var = []
#        ing = 0
#        for key, value in element.iteritems():
#            keys.append(key)
#            ing += 1
#            var.append(str(value))
#        cursor.execute(query+", ".join(keys)+query2 + ", ".join(["?"]*ing)+")", var)
#connection.commit()

from xml.dom import minidom
doc = minidom.parse("/home/jjgomera/pychemqt/images/equipment/ciclon.svg")

for input in doc.getElementsByTagName("ins")[0].childNodes:
    if isinstance(input, minidom.Element):
        if input.tagName == "in":
            x=input.getAttribute("x")
            y=input.getAttribute("y")
            d=input.getAttribute("d")
            print(dir(x))
            print(x, y, d)
            
for input in doc.getElementsByTagName("outs")[0].childNodes:
    if isinstance(input, minidom.Element):
        if input.tagName == "out":
            x=input.getAttribute("x")
            y=input.getAttribute("y")
            d=input.getAttribute("d")
            print(x, y, d)
            
#print(dir(doc.getElementsByTagName('ins')[0]))
#while input.tagName!="ins":
#    node=node.nextSibling()
#child=input.firstChild
#while child:
#    x=child.getAttribute("x")
#    print(x)
#    y=child.attribute("y")
#    d=child.attribute("d")
#    input.append([x, y, d])
#    child = child.nextSibling()
#doc.unlink()

