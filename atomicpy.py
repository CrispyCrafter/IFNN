import pandas as pd
from bs4 import BeautifulSoup
import requests
import os
import numpy as np
import json as jsonconvert



import numpy as np
import re

import random
random.random()

class AtomInfo(object):
    """Simple Class to retrieve atomic radii infromation from wikipedia """

    def __init__(self):
        dir_path = os.getcwd()
        VdW_PATH = {"path": os.path.join(dir_path,"VdW")}

        self.outfile = os.path.join(VdW_PATH["path"],"element_radii.csv")

        self.frame = self.retrieve()
        self.json = self.json()

    def json(self,query='van der Waals'):
        elements = self.frame[['symbol',query]]
        elements = elements.dropna()
        columns = elements.T.iloc[0].values
        values = elements.T.iloc[1].values

        jsonFrame = pd.DataFrame(values).T
        jsonFrame.columns = columns
        jsonString = jsonFrame.to_json(orient='records')
        return jsonconvert.loads(jsonString)[0]

    def CleanFrame(self,frame):
        frame = frame.drop(columns=[frame.columns[0]])
        frame['Metallic'].replace(regex=True,inplace=True,to_replace=r'\D',value=r'')
        frame["empirical †"].replace(regex=True,inplace=True,to_replace=r'\D',value=r'')
        frame['Calculated'].replace(regex=True,inplace=True,to_replace=r'\D',value=r'')
        frame["van der Waals"].replace(regex=True,inplace=True,to_replace=r'\D',value=r'')
        frame['Covalent (single bond)'].replace(regex=True,inplace=True,to_replace=r'\D',value=r'')
        frame["Covalent (triple bond)"].replace(regex=True,inplace=True,to_replace=r'\D',value=r'')

        frame.replace("",np.NaN)
        frame.replace("",np.NaN)

        convert_cols = ['Metallic',
                        'empirical †',
                        'Calculated',
                        'van der Waals',
                        'Covalent (single bond)',
                        'Covalent (triple bond)']

        for col in convert_cols:
            frame[col] = pd.to_numeric(frame[col],errors='coerce')


        frame[convert_cols] = frame[convert_cols]*0.01

        return frame


    def retrieve(self):
        outfile = self.outfile

        if os.path.isfile(outfile):
            frame = pd.read_csv(outfile)

            clean_frame = self.CleanFrame(frame)

            return clean_frame
        else:
            url = "https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)"

            r = requests.get(url)
            r.content

            soup = BeautifulSoup(r.content, 'html5lib')

            table = soup.find("table")
            headers = table.find_all("th")
            index = []
            for header in headers:
                index.append(header.get_text())

            rows = table.find_all("tr")
            rows.pop(0)


            elements = []
            for row in rows:
                element_row = row.find_all("td")
                element = []
                for column in element_row:
                    element.append(column.get_text())

                elements.append(element)

            element_frame = pd.DataFrame(elements)
            element_frame.columns = index

            element_frame.to_csv(path_or_buf=outfile)

            return self.CleanFrame(element_frame)

class Molecule(object):
    def __init__(self,read=None):
        if os.path.isfile(read):
            self.filename = read
            self.readfile()
            self.json_gen()

    def readfile(self):
        with open(self.filename, "r") as f:
            xyz = f.read()
            xyz = xyz.split("\n")

        digits = re.compile("\W\d+\.\d+")
        letters = re.compile(r"\b[A-Z]\b")

        self.coords = []
        self.atoms = []

        for i in xyz:
            d = digits.findall(i)
            a = letters.findall(i)
            if d:
                self.coords.append(np.array(d,np.float64))
                self.atoms.append(a[0])

    def json_gen(self):
        try:
            self.json = {"atoms": np.array(self.atoms), "coords" : np.array(self.coords)}
        except Exception as e:
            print(e)

    def write_dummy(self,Random=True):
        with open(self.filename.split(".xyz")[0] + "_trans.xyz", "w") as f:
            if Random:
                rand = np.array(random.sample(range(3,7),  1))[0]
            else:
                rand = 10
            for i, var in enumerate(self.atoms):
                if i == 0:
                    f.write("{}\n\n".format(len(self.atoms)))
                f.write('{} {:13.12f} {:13.12f} {:13.12f}\n'.format(
                    self.atoms[i][0],
                    self.coords[i][0] + rand,
                    self.coords[i][1] + rand,
                    self.coords[i][2] + rand)
                        )

    def write(self,name,newcoords):
        with open(self.filename.split(".xyz")[0] + name + ".xyz", "w") as f:
            for i, var in enumerate(self.atoms):
                if i == 0:
                    f.write("{}\n\n".format(len(self.atoms)))
                f.write('{} {:13.12f} {:13.12f} {:13.12f}\n'.format(
                    self.atoms[i][0],
                    newcoords[i][0],
                    newcoords[i][1],
                    newcoords[i][2])
                        )

    def find_centroid(self):
        self.centroid = []
        for i,val in enumerate(self.coords):
            if i == 0:
                self.centroid = self.coords[0]/(i+1)
            else:
                self.centroid = (self.centroid + self.coords[i])/(i+1)

        return np.array(self.centroid)

    def transf_centroid(self):
        self.centroid = self.find_centroid()

        self.centroid_xyz = [vec - self.centroid for vec in self.coords]

        with open(self.filename.split(".xyz")[0] + "_transf_centroid.xyz", "w") as f:
            for i, var in enumerate(self.atoms):
                if i == 0:
                    f.write("{}\n\n".format(len(self.atoms)))
                f.write('{} {:13.12f} {:13.12f} {:13.12f}\n'.format(self.atoms[i][0], self.centroid_xyz[i][0], self.centroid_xyz[i][1], self.centroid_xyz[i][2]))
        return print("Coordinates tranformed to centroid {}".format(self.centroid))
