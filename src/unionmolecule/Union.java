/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package unionmolecule;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.Bond;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.interfaces.Algorithm;
import org.openscience.smsd.tools.MoleculeSanityCheck;

/**
 *
 * @author Asad <asad@ebi.ac.uk>
 */
public class Union {

    /**
     * Return Union of two molecules
     * @param mol1 
     * @param mol2 
     * @return
     * @throws IOException
     * @throws CDKException
     */
    public List<String> join(IAtomContainer mol1, IAtomContainer mol2) throws IOException, CDKException {
        /*
         *Generate molecules from SMILES
         */
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
//        IMolecule mol1 = sp.parseSmiles("OOC1=CC=CC=C1");
//        IMolecule mol2 = sp.parseSmiles("c1ccc(cc1)c2ccccc2");

        /*
         *Set atom ids in query
         */
        int i = 0;
        for (IAtom atom1 : mol1.atoms()) {
            atom1.setID(String.valueOf((i++)));
        }

        /*
         *Set atom ids in target
         */
        int j = 0;
        for (IAtom atom2 : mol2.atoms()) {
            atom2.setID(String.valueOf((j++)));
        }

        /*
         * Configure query molecule, add missing hydrogens and aromatise 
         * rings in the molecules 
         */
        MoleculeSanityCheck.aromatizeMolecule(mol1);

        /*
         * Configure target molecule, add missing hydrogens and aromatise 
         * rings in the molecules 
         */
        MoleculeSanityCheck.aromatizeMolecule(mol2);

        /*
         * Perform atom-atom mapping between query and target.
         * Extract maximum common subgraphs.
         */
        Isomorphism isomorphism = new Isomorphism(Algorithm.DEFAULT, true);
        isomorphism.init(mol1, mol2);
        isomorphism.setChemFilters(false, false, false);

        /*
         * Solution counter
         * 
         */
        int combinations = 1;

        /*
         * This list holds SMILES of unique solutions
         * 
         */
        List<String> acSet = new ArrayList<String>();

        /*
         * Check if an isomorphism exist between query and target
         * 
         */
        if (isomorphism.getFirstAtomMapping() != null) {

            /*
             * Iterate over MCS solutions
             * 
             */
            for (Map<IAtom, IAtom> map : isomorphism.getAllAtomMapping()) {

                /* 
                 * Mathematical aim
                 * (A ∪ B) = n(A) + n(B) - (A ∩ B)
                 * Therefore, link target molecule to the query molecule.
                 * MCS/intersection between query and target is
                 * skipped and remaining target structure is added.
                 *
                 * Union Atom Container initialised
                 * 
                 */
                IAtomContainer union = new AtomContainer();

                /*
                 * Add query atoms to the union atom container
                 * 
                 */
                for (IAtom atom : mol1.atoms()) {
                    union.addAtom(atom);
                }

                /*
                 * Add query bonds to the union atom container
                 * 
                 */
                for (IBond bond : mol1.bonds()) {
                    union.addBond(bond);
                }

                /*
                 * Add target atoms and bonds to the union atom container
                 * 
                 */
                for (IBond bond : mol2.bonds()) {
                    IAtom a1 = bond.getAtom(0);
                    IAtom a2 = bond.getAtom(1);

                    /*
                     * Add target atoms and bonds which are not common
                     * 
                     */
                    if (!map.containsValue(a1)
                            && !map.containsValue(a2)) {
                        if (!union.contains(a1)) {
                            union.addAtom(a1);
                        }
                        if (!union.contains(a2)) {
                            union.addAtom(a2);
                        }
                        union.addBond(bond);
                    } else if (map.containsValue(a1)
                            && !map.containsValue(a2)) {
                        if (!union.contains(a2)) {
                            union.addAtom(a2);
                        }
                        union.addBond(new Bond(a2, getKey(a1, map), bond.getOrder(), bond.getStereo()));
                    } else if (!map.containsValue(a1)
                            && map.containsValue(a2)) {
                        if (!union.contains(a1)) {
                            union.addAtom(a1);
                        }
                        union.addBond(new Bond(a1, getKey(a2, map), bond.getOrder(), bond.getStereo()));
                    }
                }

                /*
                 * check if this combination is chemically valid
                 */
                if (isChemicallyValid(union)) {
                    String molSMILES = getSMILES(union).toString();
                    if (!acSet.contains(molSMILES)) {
                        acSet.add(molSMILES);
                    }
                }
            }
        }
//        /*
//         * Print unique union atom container solutions
//         * 
//         */
//        for (String container : acSet) {
//            System.out.println("\n-------------" + " Combination " + combinations++ + "--------------------");
//            System.out.println("Query SMILES " + getSMILES(mol1).toString() + ", count " + mol1.getAtomCount());
//            System.out.println("Target SMILES " + getSMILES(mol2).toString() + ", count " + mol2.getAtomCount());
//            System.out.println("Union SMILES " + container + ", count " + sp.parseSmiles(container).getAtomCount());
//        }
        return acSet;

    }

    /*
     * Returns mapped query atom for a give target atom
     * present in an atom-atom mapping map.
     * NULL if its not present
     *
     */
    private IAtom getKey(IAtom a1, Map<IAtom, IAtom> map) {
        for (Map.Entry<IAtom, IAtom> v : map.entrySet()) {
            if (v.getValue() == a1) {
                return v.getKey();
            }
        }
        return null;
    }

    /*
     * Returns SMILES for a molecule 
     *
     */
    /**
     * 
     * @param molecule
     * @return
     * @throws CDKException
     */
    protected String getSMILES(IAtomContainer molecule) throws CDKException {

        String smiles = "";
        if (molecule.getAtomCount() == 0) {
            return smiles;
        }

        SmilesGenerator sg = new SmilesGenerator(true);
        AllRingsFinder arf = new AllRingsFinder();
        arf.setTimeout(900000);

        IRingSet findAllRings = arf.findAllRings(molecule);
        sg.setRings(findAllRings);

        sg.setRingFinder(arf);
        smiles = sg.createSMILES(molecule, false, new boolean[molecule.getBondCount()]);
        return smiles;
    }

    /*
     * Quick and dirty check for molecule over saturation 
     * Return true is molecule is not saturated else false
     */
    /**
     * 
     * @param union
     * @return
     * @throws CDKException
     */
    protected boolean isChemicallyValid(IAtomContainer union) throws CDKException {
        for (IAtom atom : union.atoms()) {
            if ((union.getConnectedBondsCount(atom) + atom.getFormalCharge())
                    > atom.getFormalNeighbourCount()) {
                return false;
            }
        }
        return true;
    }
}