
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.time.Duration;
import java.time.Instant;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Scanner;

import org.chocosolver.solver.ResolutionPolicy;
import org.chocosolver.solver.Solver;
import org.chocosolver.solver.constraints.ICF;
import org.chocosolver.solver.constraints.IntConstraintFactory;
import org.chocosolver.solver.constraints.extension.Tuples;
import org.chocosolver.solver.propagation.hardcoded.SevenQueuesPropagatorEngine;
import org.chocosolver.solver.propagation.hardcoded.TwoBucketPropagationEngine;
import org.chocosolver.solver.search.limits.FailCounter;
import org.chocosolver.solver.search.loop.lns.LNSFactory;
import org.chocosolver.solver.search.loop.monitors.IMonitorSolution;
import org.chocosolver.solver.search.measure.IMeasures;
import org.chocosolver.solver.search.solution.BestSolutionsRecorder;
import org.chocosolver.solver.search.solution.Solution;
import org.chocosolver.solver.search.strategy.ISF;
import org.chocosolver.solver.search.strategy.IntStrategyFactory;
import org.chocosolver.solver.trace.Chatterbox;
import org.chocosolver.solver.variables.IntVar;
import org.chocosolver.solver.variables.VF;
import org.chocosolver.solver.variables.VariableFactory;
import org.chocosolver.util.ESat;

public class CryptoMain {

	private static int NBROUNDS=0;	        // The number of rounds to attack
	private static String inputFile = null;
	private static String outputFile = null;
	public int[][][] DXS1;
	public int[][][] DYS1;
	public int[][][] DKS1;
	public int[] Pk;
	public int[] Ps;
	/*
	 * AES-128:
	 * 3 rounds: 5
	 * 4 rounds: 12
	 * 5 rounds: 17
	 * 
	 * (cf Fouque)
	 */
	private int KEY_BITS=128;        // The size of the key	
	private int BLOCK_BITS=128;	    // The size of the blocks

	private int BC=BLOCK_BITS/32;    
	private int KC=KEY_BITS/32;

	static PrintWriter writer;


	public static void main(String[] args) throws FileNotFoundException {
	    inputFile=args[0];
	    outputFile=args[1];
	    NBROUNDS= Integer.parseInt(args[2]);
	    CryptoMain CM=new CryptoMain();
	    CM.Step2();
	}
	public  void Step2() {
		// Original SR and Pk
		Ps = new int[]{0,13,10,7,4,1,14,11,8,5,2,15,12,9,6,3};
		Pk = new int[]{8,1,7,15,10,4,2,3,6,9,11,0,5,12,14,13};
		
		// KLPS
		// Ps = new int[]{0,13,10,7,4,1,14,11,8,5,2,15,12,9,6,3};
		// Pk = new int[]{5,2,3,8,9,6,7,12,13,10,11,0,1,14,15,4};
		
		// perm5r
		// Ps = new int[]{0,13,10,7,4,1,14,11,8,5,2,15,12,9,6,3};
		// Pk = new int[]{15,0,2,3,4,11,5,7,6,12,8,10,9,1,13,14};

		// PsPk1
		// Ps = new int[]{0,1,2,4,3,8,9,12,5,13,14,15,6,7,10,11};
		// Pk = new int[]{10,4,12,11,6,2,5,1,8,0,9,7,13,14,15,3};

		// PsPk2
		// Ps = new int[]{0,1,2,4,3,8,9,12,5,6,13,14,7,10,11,15};
		// Pk = new int[]{15,14,11,10,6,12,4,0,3,8,1,9,2,5,13,7};

		// PsPk3
		// Ps = new int[]{0,1,4,8,9,10,12,13,5,6,14,15,2,3,7,11};
		// Pk = new int[]{14,12,8,6,7,4,0,1,3,11,10,2,9,5,13,15};

		// PsPk4
		// Ps = new int[]{0,1,2,8,4,9,12,13,5,6,7,14,3,10,11,15};
		// Pk = new int[]{12,14,11,4,8,0,3,7,10,15,2,9,6,13,5,1};

		// PsPk5
		// Ps = new int[]{0,1,2,8,4,9,12,13,3,5,14,15,6,7,10,11};
		// Pk = new int[]{5,9,15,13,3,4,6,2,11,7,10,0,8,14,1,12};

		//PsPk6
		// Ps = new int[]{0,1,2,8,4,9,12,13,5,6,7,14,3,10,11,15};
		// Pk = new int[]{14,12,11,4,8,0,3,7,10,15,2,9,6,13,5,1};

		//PsPk7
		// Ps = new int[]{0,1,2,8,4,9,12,13,5,6,7,14,3,10,11,15};
		// Pk = new int[]{12,14,11,4,10,0,3,7,8,15,2,9,6,13,5,1};

		//PsPk8
		// Ps = new int[]{0,1,2,8,4,9,12,13,5,6,7,14,3,10,11,15};
		// Pk = new int[]{12,14,11,4,8,0,3,7,2,15,10,9,6,13,5,1};

	    
		int cpt = 0,cpt2=0;
		int n=NBROUNDS;
		BufferedReader in = null;
		PrintWriter writer = null;
		String line = null;
		Scanner scanner = null;
		DXS1=new int[n][4][BC];DKS1=new int[n][4][BC];DYS1=new int[n][4][BC];
		try
		{
			in = new BufferedReader(new FileReader(inputFile));
			writer = new PrintWriter(outputFile, "UTF-8");
			scanner=new Scanner(in);
		}
		catch(FileNotFoundException | UnsupportedEncodingException exc)
		{
			exc.printStackTrace();
		}

		try {

	while ((line=in.readLine() )!= null) {
				if (line.length()>0 && line.startsWith("delta")) {
				    //System.out.println(line);
					cpt2=0;		
					String[] tab=line.substring(line.indexOf("[")+1,line.indexOf("]")).split(",");
					if (line.startsWith("deltaY")) {
						for (int r=0;r<n-1;r++) {
							for (int j=0;j<4;j++) {
								for (int i=0;i<4;i++) {
									DYS1[r][i][j]=Integer.parseInt(tab[cpt2++].trim());
								}
							}
						}
						cpt++;
					}
					else if (line.startsWith("deltaX")) {
						for (int r=0;r<n;r++) {
							for (int j=0;j<4;j++) {
								for (int i=0;i<4;i++) {
									DXS1[r][i][j]=Integer.parseInt(tab[cpt2++].trim());
								}
							}
						}
						cpt++;
					}
					else if (line.startsWith("deltaK")) {
						for (int j=0;j<4;j++) {
							for (int i=0;i<4;i++) {
								DKS1[0][i][j]=Integer.parseInt(tab[cpt2++].trim());
							}
						}
						cpt++;
						for(int r = 1; r < n ; r++){
							for(int i = 0; i < 4; i++){
								for(int j = 0; j < 4; j++){
									int pij = Pk[4*j+i];
									int j2 = pij/4;
									int i2 = pij%4;
									DKS1[r][i2][j2] = DKS1[r-1][i][j];
								}
							}
						}
						// for (int r=0;r<n;r++) {
						// 	for (int j=0;j<4;j++) {
						// 		for (int i=0;i<4;i++) {
						// 			DKS1[r][i][j]=Integer.parseInt(tab[cpt2++].trim());
						// 		}
						// 	}
						// }
						// cpt++;
					}
				}
				if (cpt==3) {
					new Step2(writer, DXS1, DYS1, DKS1, n, BC, KC, Pk, Ps);
					cpt=0;
				}
			}


		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}


