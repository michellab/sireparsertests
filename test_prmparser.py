# lohedges and @jmichel80 - Could you get example files and get sander to
# give you Amber-calculated energies so that we can double-check that the
# AmberPrm parser is reading/writing input files that have equivalent energies as
# sander. Tests could include "calculate energy in sander, then read/write the file
# using Sire, then calculate energy in sander using new file to ensure that the energy
# hasn't changed". Also, the bond, angle and dihedral+improper energies should be equal
# (dihedral+improper are summed together). The intramolecular coulomb and LJ energy can be
# compared if you use a long cutoff.
# The test_nrg function in test_amber2.py shows how to calculate intramolecular energies in
# Sire that match those calculated in Amber.

from Sire.IO import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.System import *
from Sire.Move import *
from Sire.MM import *
from Sire.FF import *
from Sire.Units import *
from Sire.Vol import *
from Sire.Maths import *
from Sire.Base import *
from Sire.Qt import *

import Sire.Stream

import glob

import os,re,sys
import shutil

from nose.tools import assert_equal, assert_almost_equal

combining_rules = "arithmetic"
temperature = 25 * celsius
pressure = 1 * atm
cutoff_dist = 1000 * angstrom

def _createSystem(molecules):
    moleculeNumbers = molecules.molNums()
    moleculeList = []

    for moleculeNumber in moleculeNumbers:
        molecule = molecules.molecule(moleculeNumber)[0].molecule()
        moleculeList.append(molecule)

    molecules = MoleculeGroup("molecules")
    ions = MoleculeGroup("ions")

    for molecule in moleculeList:
        natoms = molecule.nAtoms()
        if natoms == 1:
            ions.add(molecule)
        else:
            molecules.add(molecule)

    all = MoleculeGroup("all")
    all.add(molecules)
    all.add(ions)

    # Add these groups to the System
    system = System()

    system.add(all)
    system.add(molecules)
    system.add(ions)

    return system


def _setupForcefields(system):


    space = Cartesian()
    
    all = system[MGName("all")]
    molecules = system[MGName("molecules")]
    ions = system[MGName("ions")]

    # - first solvent-solvent coulomb/LJ (CLJ) energy
    internonbondedff = InterCLJFF("molecules:molecules")
    internonbondedff.add(molecules)

    inter_ions_nonbondedff = InterCLJFF("ions:ions")
    inter_ions_nonbondedff.add(ions)

    inter_ions_molecules_nonbondedff = InterGroupCLJFF("ions:molecules")

    inter_ions_molecules_nonbondedff.add(ions, MGIdx(0))
    inter_ions_molecules_nonbondedff.add(molecules, MGIdx(1))

    # Now solute bond, angle, dihedral energy
    intrabondedff = InternalFF("molecules-intrabonded")
    intrabondedff.add(molecules)

    # Now solute intramolecular CLJ energy
    intranonbondedff = IntraCLJFF("molecules-intranonbonded")
    intranonbondedff.add(molecules)

    # Here is the list of all forcefields
    forcefields = [internonbondedff, intrabondedff, intranonbondedff, inter_ions_nonbondedff,
                   inter_ions_molecules_nonbondedff]

    for forcefield in forcefields:
        system.add(forcefield)

    system.setProperty("space", space)
    system.setProperty("switchingFunction", CHARMMSwitchingFunction(cutoff_dist))
    system.setProperty("combiningRules", VariantProperty(combining_rules))

    total_nrg = internonbondedff.components().total() + \
                intranonbondedff.components().total() + intrabondedff.components().total() + \
                inter_ions_nonbondedff.components().total() + inter_ions_molecules_nonbondedff.components().total()
                

    e_total = system.totalComponent()

    system.setComponent(e_total, total_nrg)

    # Add a monitor that calculates the average total energy and average energy
    # deltas - we will collect both a mean average and an zwanzig average
    system.add("total_energy", MonitorComponent(e_total, Average()))

    return system

def _sanderEnergy(prm_file, rst_file):
    #sanderbin = '/home/julien/software/amber18/bin/sander'
    sanderbin = '/home/julien/local/AMBER16/amber16/bin/sander'
    protocol = """
  Single point energy calc
&cntrl
ntf = 1, 
ntb = 0, 
ntc = 1,
extdiel = 78.5,
cut = 99999.0, 
nsnb = 99999,
imin = 1, 
maxcyc = 1, 
ncyc = 0
&end
    """
    infile = open('sp.in','w')
    infile.write(protocol)
    infile.close()
    cmd = "%s -O -i sp.in -o sp.out -p %s -c %s " % (sanderbin,prm_file,rst_file)
    os.system(cmd)
    sander_out = open('sp.out','r')
    buffer = sander_out.readlines()
    bond_nrg = 0.0
    angle_nrg = 0.0
    dihed_nrg = 0.0
    vdw_nrg = 0.0
    vdw14_nrg = 0.0
    eel_nrg = 0.0
    eel14_nrg = 0.0
    for line in buffer:
        if line.find(" BOND    =") > -1:
            elems = line.split()
            bond_nrg = float(elems[2])
            angle_nrg = float(elems[5])
            dihed_nrg = float(elems[8])
        if line.find("VDWAALS =") > -1:
            elems = line.split()
            vdw_nrg = float(elems[2])
            eel_nrg = float(elems[5])
        if line.find("1-4 VDW =") > -1:
            elems = line.split()
            vdw14_nrg = float(elems[3])
            eel14_nrg = float(elems[7])
    # FUDGE EEL ENERGIES BECAUSE SANDER AND SIRE HAVE DIFFERENT VACUUM PERMITTIVITY CONSTANTS
    eel_fudge = (18.222615317)**2/(18.2223)**2
    eel_nrg = eel_fudge*eel_nrg
    eel14_nrg = eel_fudge*eel14_nrg
    tot_nrg = bond_nrg + angle_nrg+dihed_nrg+vdw_nrg+vdw14_nrg+eel_nrg+eel14_nrg
    # Tidy up
    #cmd = "rm -f sp.out sp.in restrt mdinfo"
    #os.system(cmd)
    
    return [tot_nrg*kcal_per_mol, bond_nrg, angle_nrg, dihed_nrg, vdw_nrg, vdw14_nrg, eel_nrg, eel14_nrg]

    
def test_nrg(prm_file, rst_file, ostream, verbose=False):

    if verbose:
        ostream.write("############### Testing %s %s \n"%  (prm_file,rst_file))
    # Step 1 read input in Sire. Compute sp energy.
    #system = MoleculeParser.read(prm_file, rst_file)
    #molecules = system.molecules()
    molecules, space = Amber().readCrdTop(rst_file, prm_file)
    system = _createSystem(molecules)
    system = _setupForcefields(system)
    nrg1 = system.energy()
    # Step 2. Compute sp energy with sander.
    sander_energies1 = _sanderEnergy(prm_file, rst_file)
    diff1 = nrg1-sander_energies1[0]
    #import pdb; pdb.set_trace()
    if verbose:
        ostream.write("Difference in energies between sander and sire with input %s \n" % diff1)
    assert_almost_equal( nrg1.value(), sander_energies1[0].value(), 2)
    # Tidy up sabder output
    cmd = "rm -f sp.out sp.in restrt mdinfo"
    os.system(cmd)
    # Step 3. Write input from Sire. Read again. Compute sp energy.
    # Write back to file.
    p = AmberRst7(system)
    p.writeToFile("test.rst7")
    p = AmberPrm(system)
    p.writeToFile("test.prm7")
    s2 = MoleculeParser.read(["test.rst7", "test.prm7"])
    molecules2 = s2.molecules()
    s2 = _createSystem(molecules2)
    s2 = _setupForcefields(s2)
    nrg2 = s2.energy()
    diff2 = nrg1 - nrg2
    if verbose:
        ostream.write("Difference in energies between sire input and sire output %s \n" % diff2)
    assert_almost_equal( nrg1.value(), nrg2.value(), 2)
    # Step 4. Read again written input in sander. Compute sp energy with sander.
    sander_energies2 = _sanderEnergy("test.prm7","test.rst7")
    diff3 = nrg2-sander_energies2[0]
    if verbose:
        ostream.write("Difference in energies between sire output and sander output %s \n" % diff3)
    assert_almost_equal( nrg2.value(), sander_energies2[0].value(), 2)
    # Tidy up
    cmd = "rm -f sp.out sp.in restrt mdinfo"
    os.system(cmd)
    #tidy up
    cmd = "rm test.rst7 test.prm7"
    os.system(cmd)
    #import pdb; pdb.set_trace()
    #assert_almost_equal( e_bond, 5.1844, 2 )

if __name__ == '__main__':
    #rst_files = glob.glob('input/*/03-fesetup-morph/*/*/*equili*')
    #prm_files = glob.glob('input/*/03-fesetup-morph/*/*/solvated.parm7')

    rst_files = glob.glob('input/*/03-fesetup-morph/_ligands/*/*equili*')
    prm_files = glob.glob('input/*/03-fesetup-morph/_ligands/*/solvated.parm7')

    #rst_files = glob.glob('input/*/03-fesetup-morph/_complexes/*/*equili*')
    #prm_files = glob.glob('input/*/03-fesetup-morph/_complexes/*/solvated.parm7')

    #prm_files = ['input/bace_ds/03-fesetup-morph/_ligands/bace_lig19/solvated.parm7']
    #rst_files = ['input/bace_ds/03-fesetup-morph/_ligands/bace_lig19/bace_lig19_equilibrated.rst7']

    #prm_files = ['input/bace_ds/03-fesetup-morph/_complexes/BACE:bace_lig1/solvated.parm7']
    #rst_files = ['input/bace_ds/03-fesetup-morph/_complexes/BACE:bace_lig1/BACE:bace_lig1_equilibrated.rst7']
    
    #prm_files = ['input/tyk2_ds/03-fesetup-morph/_ligands/tyk_lig15/solvated.parm7']
    #rst_files = ['input/tyk2_ds/03-fesetup-morph/_ligands/tyk_lig15/tyk_lig15_equilibrated.rst7']

    prm_files = ['input/jnk1_ds/03-fesetup-morph/_ligands/jnk1_lig17/solvated.parm7']
    rst_files = ['input/jnk1_ds/03-fesetup-morph/_ligands/jnk1_lig17/jnk1_lig17_equilibrated.rst7']

    #prm_files = ['input/tyk2_ds/03-fesetup-morph/_complexes/TYK2:tyk_lig1/solvated.parm7']
    #rst_files = ['input/tyk2_ds/03-fesetup-morph/_complexes/TYK2:tyk_lig1/TYK2:tyk_lig1_equilibrated.rst7']

    logfile = 'results.log'
    #logfile = sys.stdout
    #ostream = open(logfile,'w')
    ostream = sys.stdout

    ostream.write("Size of test set: %s \n" % len(rst_files))
    #    test_nrg(True)
    for x in range(0,len(prm_files)):
        prm = prm_files[x]
        rst = rst_files[x]
        try:
            test_nrg(prm, rst, ostream,verbose=True)
        except AssertionError:
            ostream.write("TEST FAILED! FOR %s %s \n" % (prm,rst))
        ostream.flush()
