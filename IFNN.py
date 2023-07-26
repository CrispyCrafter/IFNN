import numpy as np
from numpy import sin, cos
import re

from scipy.optimize import minimize
from scipy.optimize import basinhopping
from scipy.spatial.distance import euclidean as euclid
import random
import os


""" Importing structure still messed up with ipython kernel """
try:
        import atomicpy
except Exception as e:
        print(e)


class IFNN(object):
        """IFNN pre-optimisation class to generate starting cooridnates for ab-initio modelling """
        def __init__(self, molAclass,molBclass,VdW):
                # Initialise file dictionaries
                self.molA = molAclass.json
                self.molB = molBclass.json

                # recall molclasses()

                self.molAclass = molAclass
                self.molBclass = molBclass

                # Initialise file directories
                self.pathA = os.path.dirname(os.path.realpath(molAclass.filename))
                self.pathB = os.path.dirname(os.path.realpath(molBclass.filename))

                filenameA = os.path.basename(molAclass.filename).split(".xyz")[0]
                filenameB = os.path.basename(molBclass.filename).split(".xyz")[0]

                self.filename = os.path.join(self.pathA,filenameA+"_"+filenameB)

                self.OPT = self.genRandOpt()
                self.VdW = VdW

                self.results = np.array([])

        def genRandOpt(self):
                 OPT = np.array([
                     random.sample(range(3,10),  1)[0],
                     random.sample(range(3,10),  1)[0],
                     random.sample(range(3,10),  1)[0],
                     random.sample(range(0,360),  1)[0],
                     random.sample(range(0,360),  1)[0],
                     random.sample(range(0,360),  1)[0],
                     ])
                 return OPT

        def genRandExclOpt(self):

                excA = self.molAclass.exclusion()
                excB = self.molBclass.exclusion()
                AB_Excl = excA + excB

                polAngles = np.deg2rad([np.random.randint(0,360) for i in range(3)])

                polTrans = np.array([
                    sin(polAngles[0])*cos(polAngles[1]),
                    sin(polAngles[0])*sin(polAngles[1]),
                    cos(polAngles[0])
                    ])

                polXYZ = polTrans*AB_Excl


                return np.append(polXYZ,polAngles)

        def genTrMatrix(self,vector,OPT=None):
                if OPT is None:
                        OPT = self.OPT

                phi     = np.deg2rad(OPT[3])
                theta   = np.deg2rad(OPT[4])
                psi     = np.deg2rad(OPT[5])

                rotateZ = np.array([
                            [ cos(phi)   ,-sin(phi) , 0          ],
                            [ sin(phi)   , cos(phi) , 0          ],
                            [ 0          , 0        , 1          ]
                    ])

                rotateY = np.array([
                            [ cos(theta) , 0        , sin(theta) ],
                            [ 0          , 1        , 0          ],
                            [ -sin(theta), 0        , cos(theta) ]
                    ])

                rotateX = np.array([
                            [ 1          , 0        , 0          ],
                            [ 0          , cos(psi) ,-sin(psi)   ],
                            [ 0          , sin(psi) , cos(psi)   ]
                    ])

                translate = np.array([OPT[0], OPT[1], OPT[2]])

                rotateZYX = np.matmul(np.matmul(rotateZ,rotateY),rotateX)
                transform = np.matmul(rotateZYX,vector - translate) +  translate

                return transform

        def genConstraints(self,OPT=None):
                if OPT is None:
                        OPT = self.OPT

                VdW_cons =  np.array([[self.VdW[atomA]+ self.VdW[atomB] for atomA in self.molA["atoms"]] for atomB in self.molB["atoms"]])
                Result_cons = self.distance(OPT,False)

                nett = np.array(Result_cons - VdW_cons).flatten()
                index = nett.argsort()[:6]

                return nett[index]

        def genCoords(self,OPT=None):
                if OPT is None:
                    print("Transformation not specified.")
                    print("Initial coords of molecule A returned conditions returned")
                    OPT = self.OPT
                elif OPT:
                    OPT = self.genRandExclOpt()

                transA = np.array([ self.genTrMatrix(atom,OPT) for atom in self.molA["coords"]])
                tmp_coords = {
                        'atoms': self.molA['atoms'],
                        'coords': transA
                             }

                return tmp_coords


        def distance(self,OPT,np_sum=True,transMol=None):
                if transMol is None:
                    transMol = self.molA


                transA = np.array([ self.genTrMatrix(atom,OPT) for atom in transMol["coords"]])
                calc = [ [ euclid(atom,np.array(b)) for atom in transA ] for b in transMol["coords"] ]
                if np_sum:
                        return np.sum(calc)
                elif not np_sum:
                        return np.array(calc)
                else:
                        raise({}.format("Debugging Error"))

        def runOpt(self,OPT=None,transMol=None):
                if OPT is None:
                    OPT = self.OPT
                if transMol is None:
                    transMol = self.molA

                cons = ({'type' : 'ineq', 'fun': self.genConstraints})
                bounds = ((None,None),(None,None),(None,None),(-180,180),(None,None),(None,None))

                SLSQP = minimize(self.distance, OPT, args=(True,transMol), method='SLSQP', options={'disp': True, 'maxiter' : 1000}, constraints=cons, bounds=bounds)

                self.results = np.append(self.results,np.array([dict(SLSQP)]),axis=0)

                return dict(SLSQP)

        def runBasinOpt(self,OPT=None,transMol=None,niter=10, display=False, maxiter=1000):
                if OPT is None:
                    OPT = self.OPT
                if transMol is None:
                    transMol = self.molA

                cons = ({'type' : 'ineq', 'fun': self.genConstraints})
                bounds = ((None,None),(None,None),(None,None),(-180,180),(None,None),(None,None))

                minimizer_kwargs = {"args":(True,transMol), "method": "SLSQP", "options": {'disp': display, 'maxiter' : maxiter},'constraints':cons, 'bounds': bounds}
                basinSLSQP = basinhopping(self.distance, OPT, minimizer_kwargs=minimizer_kwargs,niter=niter)


                self.results = np.append(self.results,np.array([dict(basinSLSQP)]),axis=0)

                return dict(basinSLSQP)


        def writeXYZ(self,molecule,resultID=0):
                with open("{}_{}.xyz".format(self.filename,resultID), "w") as f:
                    for i, var in enumerate(molecule["atoms"]):
                        if i == 0:
                            f.write("{}\n\n".format(len(molecule["atoms"])))
                        f.write('{} {:13.12f} {:13.12f} {:13.12f}\n'.format(
                            molecule["atoms"][i],
                            molecule["coords"][i][0],
                            molecule["coords"][i][1],
                            molecule["coords"][i][2])
                                )

        def writeAll(self,clean=False,name=None):
                # if clean:
                #         os.
                # # Populate atom names

                dimerAtoms = np.append(self.molB['atoms'],self.molB['atoms'])
                resultNumber = 0

                print(len(dimerAtoms))

                for result in self.results:
                        if result['success']:

                            OPT = result['x']
                            dimerCoords = []
                            molA_result =  np.array([ self.genTrMatrix(atom,OPT) for atom in self.molA["coords"]])
                            dimerCoords = np.append(self.molB['coords'],molA_result, axis=0)

                            moljson = {"atoms": dimerAtoms, "coords": dimerCoords}
                            self.writeXYZ(moljson,resultNumber)
                            result_number += 1
