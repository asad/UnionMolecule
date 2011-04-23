/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package unionmolecule;

import java.io.IOException;
import java.util.List;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;

/**
 *
 * @author Asad <asad@ebi.ac.uk>
 */
public class UnionMolecule {

    /**
     * @param args the command line arguments
     * @throws InvalidSmilesException
     * @throws IOException
     * @throws CDKException  
     */
    public static void main(String[] args) throws InvalidSmilesException, IOException, CDKException {
        if (args.length < 2 || args.length > 2) {
            System.out.println("Please enter Query and Target SMILES");
            System.out.println("Example:\nOOC1=CC=CC=C1 c1ccc(cc1)c2ccccc2");
            return;
        }

        // TODO code application logic here
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer mol1 = sp.parseSmiles(args[0]);
        IAtomContainer mol2 = sp.parseSmiles(args[1]);

        List<String> acSet = new Union().join(mol1, mol2);
        int combinations = 1;
        for (String container : acSet) {
            System.out.println("\n-------------" + " Combination " + combinations++ + "--------------------");
            System.out.println("Query SMILES " + args[0].toString() + ", count " + mol1.getAtomCount());
            System.out.println("Target SMILES " + args[1].toString() + ", count " + mol2.getAtomCount());
            System.out.println("Union SMILES " + container + ", count " + sp.parseSmiles(container).getAtomCount());
        }

    }
}
