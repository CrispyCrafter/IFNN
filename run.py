""" Importing structure still messed up with ipython kernel """
try:
        import atomicpy
        import IFNN
except Exception as e:
        print(e)
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
run.opt()

print("Running 4 basin hopping  Sequential Least SQuares Programming optimisations")
run.basinopt(niter=4,display=True)

print("Writing initial and optimised result files")
run.writeAll()
