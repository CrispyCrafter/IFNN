""" Importing structure still messed up with ipython kernel """

import atomicpy
import IFNN as IFNN
import os

DIRPATH = os.getcwd()
MOLPATH = os.path.join(DIRPATH,"Molecules")

VdW = atomicpy.AtomInfo().json;

print("Reading initial xyz file")
molA = atomicpy.Molecule(os.path.join(MOLPATH,"Benzene.xyz"));

print("Creating example dummy molecule")
molA.write_dummy()
molB = atomicpy.Molecule(os.path.join(MOLPATH,"Benzene_trans.xyz"))

print("Molecules intialiased preparing optimisation parameters")

# IFFN Class takes in two 'Molecule' classes genetated from atomicpy.Molecule
# VdW is a dictionary containing the known van der Waals radii of elements

run = IFNN.IFNN(molA,molB,VdW)

print("Running single Sequential Least SQuares Programming optimisation")
x = run.runOpt()

{x}
print("Running 4 basin hopping  Sequential Least SQuares Programming optimisations")
run.runBasinOpt(niter=4,display=True)

print("Writing initial and optimised result files")
run.writeAll()

monte = [run.genCoords(True) for i in range(10)];
monte;

monteResutls = [[run.runOpt(OPT=None,transMol=monte[i])] for i in range(10)]
type(monteResutls[0][0])

minDistances = []

for result in monteResutls:
    result = dict(result[0])
    if result["success"]:
        minDistances.append([result['fun'],result['x']])

minDistances


dict(monteResutls[0][0])
run.writeAll()
run.results
