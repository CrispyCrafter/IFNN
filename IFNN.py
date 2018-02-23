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

                # Initialise file directories
                self.pathA = os.path.dirname(os.path.realpath(molAclass.filename))
                self.pathB = os.path.dirname(os.path.realpath(molBclass.filename))

                filenameA = os.path.basename(molAclass.filename).split(".xyz")[0]
                filenameB = os.path.basename(molBclass.filename).split(".xyz")[0]

                self.filename = os.path.join(self.pathA,filenameA+"_"+filenameB)

                self.OPT = self.generateOPT()
                self.VdW = VdW

                self.results = np.array([[0,0,0,0,0,0]])

        def generateOPT(self):
                 OPT = np.array([
                     random.sample(range(3,10),  1)[0],
                     random.sample(range(3,10),  1)[0],
                     random.sample(range(3,10),  1)[0],
                     random.sample(range(0,360),  1)[0],
                     random.sample(range(0,360),  1)[0],
                     random.sample(range(0,360),  1)[0],
                     ])
                 return OPT

        def generate_matrix(self,vector,OPT=None):
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

        def test(self,OPT=None):
                if OPT is None:
                        OPT = self.OPT

                VdW_test =  np.array([[self.VdW[atomA]+ self.VdW[atomB] for atomA in self.molA["atoms"]] for atomB in self.molB["atoms"]])
                Result_test = self.distance(OPT,False)

                nett = np.array(Result_test - VdW_test).flatten()
                index = nett.argsort()[:6]

                return nett[index]

        def distance(self,temp_OPT,np_sum=True):

                transA = np.array([ self.generate_matrix(atom,temp_OPT) for atom in self.molA["coords"]])
                calc = [ [ euclid(a,np.array(b)) for a in transA ] for b in self.molB["coords"] ]
                if np_sum:
                        return np.sum(calc)
                elif not np_sum:
                        return np.array(calc)
                else:
                        raise({}.format("Debugging Error"))

        def opt(self):
                cons = ({'type' : 'ineq', 'fun': self.test})
                bounds = ((None,None),(None,None),(None,None),(-180,180),(None,None),(None,None))

                SLSQP = minimize(self.distance, self.OPT, args=(True), method='SLSQP', options={'disp': True, 'maxiter' : 1000}, constraints=cons, bounds=bounds)

                self.SLS = SLSQP
                self.results = np.append(self.results,np.array([SLSQP.x]),axis=0)

                return SLSQP.x

        def basinopt(self,niter=10, display=False, maxiter=1000):
                cons = ({'type' : 'ineq', 'fun': self.test})
                bounds = ((None,None),(None,None),(None,None),(-180,180),(None,None),(None,None))

                minimizer_kwargs = {"args":(True), "method": "SLSQP", "options": {'disp': display, 'maxiter' : maxiter},'constraints':cons, 'bounds': bounds}
                basinSLSQP = basinhopping(self.distance, self.OPT, minimizer_kwargs=minimizer_kwargs,niter=niter)

                self.SLSbasin = basinSLSQP
                self.results = np.append(self.results,np.array([basinSLSQP.x]),axis=0)

                return basinSLSQP.x


        def writeXYZ(self,molecule,result_number=0):
                with open("{}_{}.xyz".format(self.filename,result_number), "w") as f:
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

                dimeratoms = np.append(self.molB['atoms'],self.molB['atoms'])
                result_number = 0

                print(len(dimeratoms))

                for result in self.results:
                        dimer_coords = []
                        molA_result =  np.array([ self.generate_matrix(atom,result) for atom in self.molA["coords"]])
                        dimer_coords = np.append(self.molB['coords'],molA_result, axis=0)

                        moljson = {"atoms": dimeratoms, "coords": dimer_coords}
                        self.writeXYZ(moljson,result_number)
                        result_number += 1
