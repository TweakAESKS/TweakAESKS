
import java.io.PrintWriter;
import org.chocosolver.samples.AbstractProblem;
import org.chocosolver.solver.ResolutionPolicy;
import org.chocosolver.solver.Solver;
import org.chocosolver.solver.constraints.ICF;
import org.chocosolver.solver.constraints.IntConstraintFactory;
import org.chocosolver.solver.constraints.LCF;
import org.chocosolver.solver.constraints.extension.Tuples;
import org.chocosolver.solver.search.solution.Solution;
import org.chocosolver.solver.search.strategy.IntStrategyFactory;
import org.chocosolver.solver.trace.Chatterbox;
import org.chocosolver.solver.variables.IntVar;
import org.chocosolver.solver.variables.VF;
import org.chocosolver.solver.variables.VariableFactory;

public class Step2 extends AbstractProblem {
	static int[][][]  shifts= {
		{  {0, 0},
			{1, 3},
			{2, 2},
			{3, 1}
		},{

			{0, 0},
			{1, 5},
			{2, 4},
			{3, 3},
		},{
			{0, 0},
			{1, 7},
			{3, 5},
			{4, 4}
		}
	}; 
	public int nbtest=0;
	public int nk;
	public int BC, KC,SC;
	public int n, best;
	public IntVar[][][] DY;
	public IntVar[][][] DX;
	public IntVar[][][] DSR;
	public IntVar[][][] DK;
	public String strategy;

	public int[] Pk;
	public int[] Ps;

	public int[][][] DYS1;
	public int[][][] DXS1;
	public int[][][] DKS1;
	public int[][] CSRYS1;

	public float totalTime;

	public IntVar[][][] DeltaSY;
	//public IntVar[][] DeltaSK;
	public PrintWriter writer;
	public IntVar[] allVar;
	public IntVar[][][] p;
	//public IntVar[][] pk;
	public IntVar[] allVarsToMax;
	public IntVar[] allProbas;

	public IntVar obj;
	public int nbSB;
	public int nbVar;	
	public double totalPXB;	

	public Step2(PrintWriter file, int[][][] X, int[][][] Y, int[][][]K, int nb, int Bc, int Kc, int[] Pk_in, int[] Ps_in){
		n = nb;

		BC=Bc;
		KC=Kc;
		nk=n*BC/KC;
		switch(BC) {
		case 4: 
			SC=0;
			break;
		case 6:
			SC=1;
			break;

		case 8: 
			SC=2;
			break;

		}

		writer= file;
		best=0;
		DYS1=Y;DKS1=K;DXS1=X;
		Pk = Pk_in; Ps = Ps_in;
		CSRYS1=new int[n][BC];
		int nS=0;

		for (int r=0;r<n;r++) {
			for (int i=0;i<4;i++) {
				for (int j=0;j<BC;j++) {
					// if ((BC*r+j)%KC==KC-1 && DKS1[r][i][j]==1) {
					// 	nS++;
					// }
					if (DXS1[r][i][j]==1) {
						nS++;
						CSRYS1[r][(j+shifts[SC][i][1])%BC]++;
					}
				}
			}
		}		        	

		totalPXB=0.0;
		nbSB=nS;
		//printStep1();
		System.out.println(nS + " Sboxes");

		// nbVar=n*4*4*BC /* X, Y, K, p */+ (n-1)*2*4*BC/* SY, SR */ + ((int)(Math.ceil((n*BC)/(double)KC))-1)*4 /* SK */+((n*BC)/KC)*4 /* pk */+1 /*obj*/;
		nbVar=n*4*4*BC /* X, Y, K, p */+ (n-1)*2*4*BC/* SY, SR */+1 /*obj*/;
		strategy = "FC";
		//strategy = "GAC2001+";
		//strategy = "AC3bit+rm";
		//strategy="STR2+";
		DY = new IntVar[n][4][BC];
		DX = new IntVar[n][4][BC];
		DK = new IntVar[n][4][BC];
		DSR = new IntVar[n-1][4][BC];
		DeltaSY= new IntVar[n-1][4][BC];
		p=new IntVar[n][4][BC];
		//DeltaSK= new IntVar[(int) (Math.ceil(n*BC/(double)KC) -1)][4];				    
		//pk=new IntVar[((n*BC)/KC)][4];

		allVar=new IntVar[nbVar];
		allVarsToMax=new IntVar[4*BC*(n)+4*((n*BC)/KC)];
		this.execute();

	}

	public void printStep1() {
		for (int j=0; j<4; j++){
			for (int k=0; k<BC; k++) System.out.print("0 ");
			System.out.print("  ");
			for (int k=0; k<BC; k++) System.out.print(DKS1[0][j][k]+" ");
			System.out.print("  ");
			for (int k=0; k<BC; k++) System.out.print(DXS1[0][j][k]+" ");
			if (j==0) {
				for (int k=0; k<BC; k++) {
					System.out.print(" "+CSRYS1[0][k]);
				}
			}
			System.out.println();
		}
		for (int i=1; i<n; i++){

			for (int j=0; j<4; j++){
				for (int k=0; k<BC; k++) System.out.print(DYS1[i-1][j][k]+" ");
				System.out.print("  ");
				for (int k=0; k<BC; k++) System.out.print(DKS1[i][j][k]+" ");
				System.out.print("  ");
				for (int k=0; k<BC; k++) System.out.print(DXS1[i][j][k]+" ");
				if (j==0) {
					for (int k=0; k<BC; k++) {
						System.out.print(" "+CSRYS1[i][k]);
					}
				}
				System.out.println();

			}

		}
		System.out.println("-------------------------------------------------");
	}



	public static int[] Sbox = new int[] {
		99, 124, 119, 123, 242, 107, 111, 197,  48,   1, 103,  43, 254, 215, 171, 118, 
		202, 130, 201, 125, 250,  89,  71, 240, 173, 212, 162, 175, 156, 164, 114, 192, 
		183, 253, 147,  38,  54,  63, 247, 204,  52, 165, 229, 241, 113, 216,  49,  21, 
		4, 199,  35, 195,  24, 150,   5, 154,   7,  18, 128, 226, 235,  39, 178, 117, 
		9, 131,  44,  26,  27, 110,  90, 160,  82,  59, 214, 179,  41, 227,  47, 132, 
		83, 209,   0, 237,  32, 252, 177,  91, 106, 203, 190,  57,  74,  76,  88, 207, 
		208, 239, 170, 251,  67,  77,  51, 133,  69, 249,   2, 127,  80,  60, 159, 168, 
		81, 163,  64, 143, 146, 157,  56, 245, 188, 182, 218,  33,  16, 255, 243, 210, 
		205,  12,  19, 236,  95, 151,  68,  23, 196, 167, 126,  61, 100,  93,  25, 115, 
		96, 129,  79, 220,  34,  42, 144, 136,  70, 238, 184,  20, 222,  94,  11, 219, 
		224,  50,  58,  10,  73,   6,  36,  92, 194, 211, 172,  98, 145, 149, 228, 121, 
		231, 200,  55, 109, 141, 213,  78, 169, 108,  86, 244, 234, 101, 122, 174,   8, 
		186, 120,  37,  46,  28, 166, 180, 198, 232, 221, 116,  31,  75, 189, 139, 138, 
		112,  62, 181, 102,  72,   3, 246,  14,  97,  53,  87, 185, 134, 193,  29, 158, 
		225, 248, 152,  17, 105, 217, 142, 148, 155,  30, 135, 233, 206,  85,  40, 223, 
		140, 161, 137,  13, 191, 230,  66, 104,  65, 153,  45,  15, 176,  84, 187,  22
	};
	public static Tuples tupleSB = Step2.createRelationSbox();
	public static Tuples createRelationSbox(){

		int trans[][] = new int[256][127];
		int probas[][] = new int[256][256];
		int ctrans[] = new int[256];
		Tuples tuples = new Tuples(true);

		for (int i=0;i<256;i++) {
			ctrans[i]=0;
			for (int j=0;j<256;j++) {
				probas[i][j]=0;
				if(j<127)
					trans[i][j]=0;
			}
		}
		for(int i=0;i<256;i++)
		{
			for(int j=0;j<256;j++)
			{
				probas[i][(Sbox[j]^Sbox[j^i])] ++;
				if(probas[i][(Sbox[j]^Sbox[j^i])]==1) {
					trans[i][ctrans[i]++]=Sbox[j]^Sbox[i^j];
				}
			}
		}
		tuples.add(0,0,0);
		int p=0;
		for (int i=1; i<256; i++){
			for (int j=0; j<127; j++) {
				p=probas[i][trans[i][j]]/2;		
				tuples.add(i,trans[i][j],p);
			//	if (i==143 && trans[i][j]==131) System.out.println("143->131 :"+ probas[i][trans[i][j]]);
			//	if (i==143 && trans[i][j]==204) System.out.println("143->204 :" + probas[i][trans[i][j]]);
			//	if (i==138 && trans[i][j]==204) System.out.println("138->204 :" + probas[i][trans[i][j]]);
			//	if (i==5 && trans[i][j]==79) System.out.println("5->79 :" + probas[i][trans[i][j]]);

			}
		}
		return tuples;
	}


	public static Tuples tupleXor = Step2.createRelationXor();
	public static Tuples createRelationXor(){
		Tuples tuples = new Tuples(true);
				for (int i=0; i<256; i++)
					for (int j=0; j<256; j++)
						tuples.add(i,j,i^j);
				return tuples;
	}
	
	public static Tuples tupleMul2xorMul3 = Step2.createRelationMul2xorMul3();
	public static Tuples createRelationMul2xorMul3(){
		// returns { (x,y,z) | x = Mul2*y xor Mul3*z }
		Tuples tuples = new Tuples(true);
		for (int i=0; i<256; i++){
			int ii;
			if (i<128) ii = 2*i;
			else ii = (2*i % 256) ^ 27;
			for (int j=0; j<256; j++){
				int jj;
				if (j<128) jj = (2*j) ^ j;
				else jj = ((2*j % 256) ^ 27) ^ j;
				tuples.add(ii ^ jj,i,j);

			}
		}
		return tuples;
	}

	public void postXorByte(IntVar x, IntVar y, IntVar z, Solver solver){// x xor y = z
		solver.post(IntConstraintFactory.table(new IntVar[]{x,y,z}, tupleXor, strategy));
	}

	

	public void postARK(IntVar[][]x, IntVar[][] y, IntVar[][] z, Solver solver){
		for (int i=0; i<4; i++){
			for (int j=0; j<BC; j++){
				postXorByte(x[i][j],y[i][j],z[i][j],solver);
			}
		}
	}


	public void postSBR0(IntVar[][] in, IntVar[][] out,IntVar[][] p, Solver solver) {


		for (int i=0; i<4; i++){
			for (int j=0;j<BC;j++) {
				LCF.ifThenElse(
						ICF.arithm(in[i][j],"=",0),
						LCF.and(
								ICF.arithm(out[i][j],"=",0),
								ICF.arithm(p[i][j],"=",0)
								),
								LCF.and(
										ICF.arithm(p[i][j],"=",2),
										ICF.table(new IntVar[]{in[i][j], out[i][j],p[i][j]}, tupleSB, strategy))
						);
			}

		}	
	}

	public void printBefore(PrintWriter writer){

		for (int j=0; j<4; j++){
			for (int k=0; k<BC; k++) writer.print(DY[0][j][k].getDomainSize() +" ");
			writer.print("  ");
			for (int k=0; k<BC; k++) writer.print(DK[0][j][k].getDomainSize()+" ");
			writer.print("  ");
			for (int k=0; k<BC; k++) writer.print(DX[0][j][k].getDomainSize()+" ");
			writer.print("  ");

			for (int k=0; k<BC; k++) writer.print(DSR[0][j][k].getDomainSize()+" ");

			writer.println();
		}
		for (int i=1; i<n; i++){
			for (int j=0; j<4; j++){
				for (int k=0; k<BC; k++) writer.print(DY[i][j][k].getDomainSize()+" ");
				writer.print("  ");
				for (int k=0; k<BC; k++) writer.print(DK[i][j][k].getDomainSize()+" ");
				writer.print("  ");
				//if (i<n-1)
				for (int k=0; k<BC; k++) writer.print(DX[i][j][k].getDomainSize()+" ");
				writer.print("  ");
				if (i<n-2)
					for (int k=0; k<BC; k++) writer.print(DSR[i][j][k].getDomainSize()+" ");
				writer.println();
			}
		}

	}








	public void postMC(IntVar[][] SR, IntVar[][] Y, Solver solver) {
		for (int i=0; i<BC; i++){
			// Y[0][i] = MUL2*SR[0][i] xor MUL3*SR[1][i] xor SR[2][i] xor SR[3][i]
			// Y[1][i] = SR[0][i] xor MUL2*SR[1][i] xor MUL3*SR[2][i] xor SR[3][i]
			// Y[2][i] = SR[0][i] xor SR[1][i] xor MUL2*SR[2][i] xor MUL3*SR[3][i]
			// Y[3][i] = MUL3*SR[0][i] xor SR[1][i] xor SR[2][i] xor MUL2*SR[3][i]
			IntVar v1 = VariableFactory.bounded("v1", 0, 255, solver);
			IntVar v2 = VariableFactory.bounded("v2", 0, 255, solver);
			IntVar v3 = VariableFactory.bounded("v3", 0, 255, solver);
			IntVar v4 = VariableFactory.bounded("v4", 0, 255, solver);
			IntVar v5 = VariableFactory.bounded("v5", 0, 255, solver);
			IntVar v6 = VariableFactory.bounded("v6", 0, 255, solver);
			IntVar v7 = VariableFactory.bounded("v7", 0, 255, solver);
			IntVar v8 = VariableFactory.bounded("v8", 0, 255, solver);	     

			solver.post(IntConstraintFactory.table(new IntVar[]{v1, SR[0][i], SR[1][i]}, Step2.tupleMul2xorMul3, strategy)); // v1 = MUL2*SR[0][i] xor MUL3*SR[1][i]
			postXorByte(v1,SR[2][i],v2,solver); // v2 = v1 xor SR[2][i]
			postXorByte(v2,SR[3][i],Y[0][i],solver); // Y[0][i] = v2 xor SR[3][i]

			solver.post(IntConstraintFactory.table(new IntVar[]{v3, SR[1][i], SR[2][i]}, Step2.tupleMul2xorMul3, strategy)); // v3 = MUL2*SR[1][i] xor MUL3*SR[2][i]
			postXorByte(v3,SR[0][i],v4,solver); // v4 = v3 xor SR[0][i]
			postXorByte(v4,SR[3][i],Y[1][i],solver); // Y[1][i] = v4 xor SR[3][i]
			solver.post(IntConstraintFactory.table(new IntVar[]{v5, SR[2][i], SR[3][i]}, Step2.tupleMul2xorMul3, strategy)); // v5 = MUL2*SR[2][i] xor MUL3*SR[3][i]
			postXorByte(v5,SR[0][i],v6,solver); // v6 = v5 xor SR[2][i]
			postXorByte(v6,SR[1][i],Y[2][i],solver); // Y[2][i] = v6 xor SR[1][i]
			solver.post(IntConstraintFactory.table(new IntVar[]{v7, SR[3][i], SR[0][i]}, Step2.tupleMul2xorMul3, strategy)); // v7 = MUL2*SR[3][i] xor MUL3*SR[0][i]
			postXorByte(v7,SR[1][i],v8,solver); // v8 = v7 xor SR[1][i]
			postXorByte(v8,SR[2][i],Y[3][i],solver); // Y[3][i] = v8 xor SR[2][i]
		}

	}

	public void postSR(IntVar[][] xIn, IntVar[][] xOut, Solver solver) {
		// for (int k=0; k<BC; k++){
		// 	solver.post(IntConstraintFactory.arithm(xIn[0][(k+shifts[SC][0][0])%BC], "=", xOut[0][k]));
		// 	solver.post(IntConstraintFactory.arithm(xIn[1][(k+shifts[SC][1][0])%BC], "=", xOut[1][k]));
		// 	solver.post(IntConstraintFactory.arithm(xIn[2][(k+shifts[SC][2][0])%BC], "=", xOut[2][k]));
		// 	solver.post(IntConstraintFactory.arithm(xIn[3][(k+shifts[SC][3][0])%BC], "=", xOut[3][k]));			
		// }
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				int pij = Ps[4*j+i];
				int j2 = pij/4;
				int i2 = pij%4;
				solver.post(IntConstraintFactory.arithm(xIn[i][j], "=", xOut[i2][j2]));
			}
		}
	}


	public void postSB(IntVar[][] xIn, IntVar[][] xOut, IntVar[][] p, Solver solver) {
		for (int i=0; i<4; i++){
			for (int j=0;j<BC;j++) {
				LCF.ifThenElse(
						ICF.arithm(xIn[i][j],"=",0),
						LCF.and(
								ICF.arithm(xOut[i][j],"=",0),
								ICF.arithm(p[i][j],"=",0)
								),
								ICF.table(new IntVar[]{xIn[i][j], xOut[i][j],p[i][j]}, tupleSB, strategy)
						);
			}

		}
	}

	public void postMCInv(IntVar[][] DX, IntVar[][] SR,  Solver solver) {

		for (int i=0; i<BC; i++){
			IntVar[] tmp=new IntVar[]{DX[0][i],DX[1][i],DX[2][i],DX[3][i]};
			IntVar[] tmp2=new IntVar[]{SR[0][i],SR[1][i],SR[2][i],SR[3][i]};

			LCF.ifThen(
					ICF.sum(tmp,"=",VF.fixed(0,solver)),
					ICF.sum(tmp2,"=",VF.fixed(0,solver))
					);


			IntVar[] tmpMC=VF.boundedArray("tmpMC", 24, 0, 255, solver);

			solver.post(IntConstraintFactory.table(new IntVar[]{DX[0][i], tmpMC[0]}, Step2.mul14, strategy)); // v1 = MUL2*SR[0][i] xor MUL3*SR[1][i]
			solver.post(IntConstraintFactory.table(new IntVar[]{DX[1][i], tmpMC[1]}, Step2.mul11, strategy)); // v1 = MUL2*SR[0][i] xor MUL3*SR[1][i]
			solver.post(IntConstraintFactory.table(new IntVar[]{DX[2][i], tmpMC[2]}, Step2.mul13, strategy)); // v1 = MUL2*SR[0][i] xor MUL3*SR[1][i]
			solver.post(IntConstraintFactory.table(new IntVar[]{DX[3][i], tmpMC[3]}, Step2.mul9, strategy)); // v1 = MUL2*SR[0][i] xor MUL3*SR[1][i]

			solver.post(IntConstraintFactory.table(new IntVar[]{DX[0][i], tmpMC[4]}, Step2.mul9, strategy)); // v1 = MUL2*SR[0][i] xor MUL3*SR[1][i]
			solver.post(IntConstraintFactory.table(new IntVar[]{DX[1][i], tmpMC[5]}, Step2.mul14, strategy)); // v1 = MUL2*SR[0][i] xor MUL3*SR[1][i]
			solver.post(IntConstraintFactory.table(new IntVar[]{DX[2][i], tmpMC[6]}, Step2.mul11, strategy)); // v1 = MUL2*SR[0][i] xor MUL3*SR[1][i]
			solver.post(IntConstraintFactory.table(new IntVar[]{DX[3][i], tmpMC[7]}, Step2.mul13, strategy)); // v1 = MUL2*SR[0][i] xor MUL3*SR[1][i]

			solver.post(IntConstraintFactory.table(new IntVar[]{DX[0][i], tmpMC[8]}, Step2.mul13, strategy)); // v1 = MUL2*SR[0][i] xor MUL3*SR[1][i]
			solver.post(IntConstraintFactory.table(new IntVar[]{DX[1][i], tmpMC[9]}, Step2.mul9, strategy)); // v1 = MUL2*SR[0][i] xor MUL3*SR[1][i]
			solver.post(IntConstraintFactory.table(new IntVar[]{DX[2][i], tmpMC[10]}, Step2.mul14, strategy)); // v1 = MUL2*SR[0][i] xor MUL3*SR[1][i]
			solver.post(IntConstraintFactory.table(new IntVar[]{DX[3][i], tmpMC[11]}, Step2.mul11, strategy)); // v1 = MUL2*SR[0][i] xor MUL3*SR[1][i]

			solver.post(IntConstraintFactory.table(new IntVar[]{DX[0][i], tmpMC[12]}, Step2.mul11, strategy)); // v1 = MUL2*SR[0][i] xor MUL3*SR[1][i]
			solver.post(IntConstraintFactory.table(new IntVar[]{DX[1][i], tmpMC[13]}, Step2.mul13, strategy)); // v1 = MUL2*SR[0][i] xor MUL3*SR[1][i]
			solver.post(IntConstraintFactory.table(new IntVar[]{DX[2][i], tmpMC[14]}, Step2.mul9, strategy)); // v1 = MUL2*SR[0][i] xor MUL3*SR[1][i]
			solver.post(IntConstraintFactory.table(new IntVar[]{DX[3][i], tmpMC[15]}, Step2.mul14, strategy)); // v1 = MUL2*SR[0][i] xor MUL3*SR[1][i]

			//solver.post(IntConstraintFactory.table(new IntVar[]{tmpMC[0], tmpMC[1]}, Step2.mul14, strategy)); // v1 = MUL2*SR[0][i] xor MUL3*SR[1][i]
			postXorByte(tmpMC[0],tmpMC[1],tmpMC[16],solver); 
			postXorByte(tmpMC[2],tmpMC[3],tmpMC[17],solver); 
			postXorByte(tmpMC[16],tmpMC[17],SR[0][i],solver); 

			postXorByte(tmpMC[4],tmpMC[5],tmpMC[18],solver); 
			postXorByte(tmpMC[6],tmpMC[7],tmpMC[19],solver); 
			postXorByte(tmpMC[18],tmpMC[19],SR[1][i],solver); 

			postXorByte(tmpMC[8],tmpMC[9],tmpMC[20],solver); 
			postXorByte(tmpMC[10],tmpMC[11],tmpMC[21],solver); 
			postXorByte(tmpMC[20],tmpMC[21],SR[2][i],solver); 

			postXorByte(tmpMC[12],tmpMC[13],tmpMC[22],solver); 
			postXorByte(tmpMC[14],tmpMC[15],tmpMC[23],solver); 
			postXorByte(tmpMC[22],tmpMC[23],SR[3][i],solver); 

		}
	}
	public static Tuples mul9=initMul9();
	public static Tuples initMul9() {
		Tuples ret=new Tuples(true);
		int tab[]= {0x00,0x09,0x12,0x1b,0x24,0x2d,0x36,0x3f,0x48,0x41,0x5a,0x53,0x6c,0x65,0x7e,0x77,
				0x90,0x99,0x82,0x8b,0xb4,0xbd,0xa6,0xaf,0xd8,0xd1,0xca,0xc3,0xfc,0xf5,0xee,0xe7,
				0x3b,0x32,0x29,0x20,0x1f,0x16,0x0d,0x04,0x73,0x7a,0x61,0x68,0x57,0x5e,0x45,0x4c,
				0xab,0xa2,0xb9,0xb0,0x8f,0x86,0x9d,0x94,0xe3,0xea,0xf1,0xf8,0xc7,0xce,0xd5,0xdc,
				0x76,0x7f,0x64,0x6d,0x52,0x5b,0x40,0x49,0x3e,0x37,0x2c,0x25,0x1a,0x13,0x08,0x01,
				0xe6,0xef,0xf4,0xfd,0xc2,0xcb,0xd0,0xd9,0xae,0xa7,0xbc,0xb5,0x8a,0x83,0x98,0x91,
				0x4d,0x44,0x5f,0x56,0x69,0x60,0x7b,0x72,0x05,0x0c,0x17,0x1e,0x21,0x28,0x33,0x3a,
				0xdd,0xd4,0xcf,0xc6,0xf9,0xf0,0xeb,0xe2,0x95,0x9c,0x87,0x8e,0xb1,0xb8,0xa3,0xaa,
				0xec,0xe5,0xfe,0xf7,0xc8,0xc1,0xda,0xd3,0xa4,0xad,0xb6,0xbf,0x80,0x89,0x92,0x9b,
				0x7c,0x75,0x6e,0x67,0x58,0x51,0x4a,0x43,0x34,0x3d,0x26,0x2f,0x10,0x19,0x02,0x0b,
				0xd7,0xde,0xc5,0xcc,0xf3,0xfa,0xe1,0xe8,0x9f,0x96,0x8d,0x84,0xbb,0xb2,0xa9,0xa0,
				0x47,0x4e,0x55,0x5c,0x63,0x6a,0x71,0x78,0x0f,0x06,0x1d,0x14,0x2b,0x22,0x39,0x30,
				0x9a,0x93,0x88,0x81,0xbe,0xb7,0xac,0xa5,0xd2,0xdb,0xc0,0xc9,0xf6,0xff,0xe4,0xed,
				0x0a,0x03,0x18,0x11,0x2e,0x27,0x3c,0x35,0x42,0x4b,0x50,0x59,0x66,0x6f,0x74,0x7d,
				0xa1,0xa8,0xb3,0xba,0x85,0x8c,0x97,0x9e,0xe9,0xe0,0xfb,0xf2,0xcd,0xc4,0xdf,0xd6,
				0x31,0x38,0x23,0x2a,0x15,0x1c,0x07,0x0e,0x79,0x70,0x6b,0x62,0x5d,0x54,0x4f,0x46};

		for (int i=0;i<256;i++) {
			ret.add(i,tab[i]);
		}

		return ret;
	}
	public static Tuples mul11=initMul11();
	public static Tuples initMul11() {
		Tuples ret=new Tuples(true);
		int tab[]= {0x00,0x0b,0x16,0x1d,0x2c,0x27,0x3a,0x31,0x58,0x53,0x4e,0x45,0x74,0x7f,0x62,0x69,
				0xb0,0xbb,0xa6,0xad,0x9c,0x97,0x8a,0x81,0xe8,0xe3,0xfe,0xf5,0xc4,0xcf,0xd2,0xd9,
				0x7b,0x70,0x6d,0x66,0x57,0x5c,0x41,0x4a,0x23,0x28,0x35,0x3e,0x0f,0x04,0x19,0x12,
				0xcb,0xc0,0xdd,0xd6,0xe7,0xec,0xf1,0xfa,0x93,0x98,0x85,0x8e,0xbf,0xb4,0xa9,0xa2,
				0xf6,0xfd,0xe0,0xeb,0xda,0xd1,0xcc,0xc7,0xae,0xa5,0xb8,0xb3,0x82,0x89,0x94,0x9f,
				0x46,0x4d,0x50,0x5b,0x6a,0x61,0x7c,0x77,0x1e,0x15,0x08,0x03,0x32,0x39,0x24,0x2f,
				0x8d,0x86,0x9b,0x90,0xa1,0xaa,0xb7,0xbc,0xd5,0xde,0xc3,0xc8,0xf9,0xf2,0xef,0xe4,
				0x3d,0x36,0x2b,0x20,0x11,0x1a,0x07,0x0c,0x65,0x6e,0x73,0x78,0x49,0x42,0x5f,0x54,
				0xf7,0xfc,0xe1,0xea,0xdb,0xd0,0xcd,0xc6,0xaf,0xa4,0xb9,0xb2,0x83,0x88,0x95,0x9e,
				0x47,0x4c,0x51,0x5a,0x6b,0x60,0x7d,0x76,0x1f,0x14,0x09,0x02,0x33,0x38,0x25,0x2e,
				0x8c,0x87,0x9a,0x91,0xa0,0xab,0xb6,0xbd,0xd4,0xdf,0xc2,0xc9,0xf8,0xf3,0xee,0xe5,
				0x3c,0x37,0x2a,0x21,0x10,0x1b,0x06,0x0d,0x64,0x6f,0x72,0x79,0x48,0x43,0x5e,0x55,
				0x01,0x0a,0x17,0x1c,0x2d,0x26,0x3b,0x30,0x59,0x52,0x4f,0x44,0x75,0x7e,0x63,0x68,
				0xb1,0xba,0xa7,0xac,0x9d,0x96,0x8b,0x80,0xe9,0xe2,0xff,0xf4,0xc5,0xce,0xd3,0xd8,
				0x7a,0x71,0x6c,0x67,0x56,0x5d,0x40,0x4b,0x22,0x29,0x34,0x3f,0x0e,0x05,0x18,0x13,
				0xca,0xc1,0xdc,0xd7,0xe6,0xed,0xf0,0xfb,0x92,0x99,0x84,0x8f,0xbe,0xb5,0xa8,0xa3};

		for (int i=0;i<256;i++) {
			ret.add(i,tab[i]);
		}

		return ret;
	}	public static Tuples mul13=initMul13();
	public static Tuples initMul13() {
		Tuples ret=new Tuples(true);
		int tab[]= {0x00,0x0d,0x1a,0x17,0x34,0x39,0x2e,0x23,0x68,0x65,0x72,0x7f,0x5c,0x51,0x46,0x4b,
				0xd0,0xdd,0xca,0xc7,0xe4,0xe9,0xfe,0xf3,0xb8,0xb5,0xa2,0xaf,0x8c,0x81,0x96,0x9b,
				0xbb,0xb6,0xa1,0xac,0x8f,0x82,0x95,0x98,0xd3,0xde,0xc9,0xc4,0xe7,0xea,0xfd,0xf0,
				0x6b,0x66,0x71,0x7c,0x5f,0x52,0x45,0x48,0x03,0x0e,0x19,0x14,0x37,0x3a,0x2d,0x20,
				0x6d,0x60,0x77,0x7a,0x59,0x54,0x43,0x4e,0x05,0x08,0x1f,0x12,0x31,0x3c,0x2b,0x26,
				0xbd,0xb0,0xa7,0xaa,0x89,0x84,0x93,0x9e,0xd5,0xd8,0xcf,0xc2,0xe1,0xec,0xfb,0xf6,
				0xd6,0xdb,0xcc,0xc1,0xe2,0xef,0xf8,0xf5,0xbe,0xb3,0xa4,0xa9,0x8a,0x87,0x90,0x9d,
				0x06,0x0b,0x1c,0x11,0x32,0x3f,0x28,0x25,0x6e,0x63,0x74,0x79,0x5a,0x57,0x40,0x4d,
				0xda,0xd7,0xc0,0xcd,0xee,0xe3,0xf4,0xf9,0xb2,0xbf,0xa8,0xa5,0x86,0x8b,0x9c,0x91,
				0x0a,0x07,0x10,0x1d,0x3e,0x33,0x24,0x29,0x62,0x6f,0x78,0x75,0x56,0x5b,0x4c,0x41,
				0x61,0x6c,0x7b,0x76,0x55,0x58,0x4f,0x42,0x09,0x04,0x13,0x1e,0x3d,0x30,0x27,0x2a,
				0xb1,0xbc,0xab,0xa6,0x85,0x88,0x9f,0x92,0xd9,0xd4,0xc3,0xce,0xed,0xe0,0xf7,0xfa,
				0xb7,0xba,0xad,0xa0,0x83,0x8e,0x99,0x94,0xdf,0xd2,0xc5,0xc8,0xeb,0xe6,0xf1,0xfc,
				0x67,0x6a,0x7d,0x70,0x53,0x5e,0x49,0x44,0x0f,0x02,0x15,0x18,0x3b,0x36,0x21,0x2c,
				0x0c,0x01,0x16,0x1b,0x38,0x35,0x22,0x2f,0x64,0x69,0x7e,0x73,0x50,0x5d,0x4a,0x47,
				0xdc,0xd1,0xc6,0xcb,0xe8,0xe5,0xf2,0xff,0xb4,0xb9,0xae,0xa3,0x80,0x8d,0x9a,0x97};

		for (int i=0;i<256;i++) {
			ret.add(i,tab[i]);
		}

		return ret;
	}	public static Tuples mul14=initMul14();
	public static Tuples initMul14() {
		Tuples ret=new Tuples(true);
		int tab[]= {0x00,0x0e,0x1c,0x12,0x38,0x36,0x24,0x2a,0x70,0x7e,0x6c,0x62,0x48,0x46,0x54,0x5a,
				0xe0,0xee,0xfc,0xf2,0xd8,0xd6,0xc4,0xca,0x90,0x9e,0x8c,0x82,0xa8,0xa6,0xb4,0xba,
				0xdb,0xd5,0xc7,0xc9,0xe3,0xed,0xff,0xf1,0xab,0xa5,0xb7,0xb9,0x93,0x9d,0x8f,0x81,
				0x3b,0x35,0x27,0x29,0x03,0x0d,0x1f,0x11,0x4b,0x45,0x57,0x59,0x73,0x7d,0x6f,0x61,
				0xad,0xa3,0xb1,0xbf,0x95,0x9b,0x89,0x87,0xdd,0xd3,0xc1,0xcf,0xe5,0xeb,0xf9,0xf7,
				0x4d,0x43,0x51,0x5f,0x75,0x7b,0x69,0x67,0x3d,0x33,0x21,0x2f,0x05,0x0b,0x19,0x17,
				0x76,0x78,0x6a,0x64,0x4e,0x40,0x52,0x5c,0x06,0x08,0x1a,0x14,0x3e,0x30,0x22,0x2c,
				0x96,0x98,0x8a,0x84,0xae,0xa0,0xb2,0xbc,0xe6,0xe8,0xfa,0xf4,0xde,0xd0,0xc2,0xcc,
				0x41,0x4f,0x5d,0x53,0x79,0x77,0x65,0x6b,0x31,0x3f,0x2d,0x23,0x09,0x07,0x15,0x1b,
				0xa1,0xaf,0xbd,0xb3,0x99,0x97,0x85,0x8b,0xd1,0xdf,0xcd,0xc3,0xe9,0xe7,0xf5,0xfb,
				0x9a,0x94,0x86,0x88,0xa2,0xac,0xbe,0xb0,0xea,0xe4,0xf6,0xf8,0xd2,0xdc,0xce,0xc0,
				0x7a,0x74,0x66,0x68,0x42,0x4c,0x5e,0x50,0x0a,0x04,0x16,0x18,0x32,0x3c,0x2e,0x20,
				0xec,0xe2,0xf0,0xfe,0xd4,0xda,0xc8,0xc6,0x9c,0x92,0x80,0x8e,0xa4,0xaa,0xb8,0xb6,
				0x0c,0x02,0x10,0x1e,0x34,0x3a,0x28,0x26,0x7c,0x72,0x60,0x6e,0x44,0x4a,0x58,0x56,
				0x37,0x39,0x2b,0x25,0x0f,0x01,0x13,0x1d,0x47,0x49,0x5b,0x55,0x7f,0x71,0x63,0x6d,
				0xd7,0xd9,0xcb,0xc5,0xef,0xe1,0xf3,0xfd,0xa7,0xa9,0xbb,0xb5,0x9f,0x91,0x83,0x8d};

		for (int i=0;i<256;i++) {
			ret.add(i,tab[i]);
		}

		return ret;
	}



	@Override
	public void buildModel() {
		long startTime = System.currentTimeMillis();

		int r,j,indP, cproba=0;
		allProbas=new IntVar[n*BC*4];

		for (int J=0;J<BC*n;J++) {
			r=J/BC;
			j=J%BC;
			indP=J/KC;
			int rk=J/KC;
			for (int i=0;i<4;i++) {
				if (r<n-1){					
					DSR[r][i][j] = VariableFactory.bounded("DeltaSR["+r+"]["+i+"]["+j+"]", 0, 255, solver);
					if (DXS1[r][i][j] == 0) DeltaSY[r][i][j] = VariableFactory.fixed(0, solver);
					else DeltaSY[r][i][j] = VariableFactory.bounded("DeltaSY["+(r)+"]["+i+"]["+j+"]", 1, 255, solver);
					if (DYS1[r][i][j] == 0){
						DY[r+1][i][j] = VariableFactory.fixed(0, solver);
					}
					else{
						DY[r+1][i][j] = VariableFactory.bounded("DeltaX["+(r+1)+"]["+i+"]["+j+"]", 1, 255, solver);
					}
				}


				if (DKS1[r][i][j] == 0){
					DK[r][i][j] = VariableFactory.fixed(0, solver);

					// if (J%KC==KC-1 && J<BC*n-1) {
					// 	DeltaSK[indP][(i+3)%4]=VF.fixed(0, solver);
					// 	pk[indP][(i+3)%4]=VF.fixed(0, solver);
					// }
					// else if (J%KC==KC-1) {
					// 	pk[indP][(i+3)%4]=VF.fixed(0, solver);
					// }
				}

				else{
					DK[r][i][j] = VariableFactory.bounded("DeltaK["+r+"]["+i+"]["+j+"]", 1, 255, solver);
					// if (J%KC==KC-1 && J<BC*n-1 && J/KC>0) {
					// 	DeltaSK[indP][(i+3)%4]=VariableFactory.bounded("DeltaSK["+r+"]["+i+"]", 1, 255, solver);
					// 	pk[indP][(i+3)%4]=VF.fixed(0, solver);
					// }
					// else if (J%KC==KC-1 && J/KC==0) {

					// 	DeltaSK[indP][(i+3)%4]=VariableFactory.bounded("DeltaSK["+r+"]["+i+"]", 1, 255, solver);
					// 	pk[indP][(i+3)%4]=VF.fixed(0, solver);
					// }
					// else if (J%KC==KC-1) 
					// 	pk[indP][(i+3)%4]=VF.fixed(0, solver);

				}
				if (DXS1[r][i][j] == 0){
					DX[r][i][j] = VariableFactory.fixed(0, solver);
					//	if (i<n-1)
					p[r][i][j]=VF.fixed(0, solver);
				}
				else {
					DX[r][i][j] = VariableFactory.bounded("DeltaY["+(r+1)+"]["+i+"]["+j+"]", 1, 255, solver);
					if (r==0)
						p[r][i][j]=VF.fixed("p["+r+"]["+i+"]["+j+"]",2, solver);
					else if (r<n-1)
						p[r][i][j]=VF.enumerated("p["+r+"]["+i+"]["+j+"]",new int[]{1,2}, solver);
					else
						p[r][i][j]=VF.fixed("p["+r+"]["+i+"]["+j+"]",2, solver);

				}			
			}
		}

		DY[0]=VF.boundedMatrix("DeltaX0", 4, BC, 0, 255,solver); 
		postSBR0(DX[0],DeltaSY[0],p[0], solver);

		for(int i=1;i<n-1;i++) {
			postSB(DX[i],DeltaSY[i],p[i], solver);
		}
		for (int i=0;i<n-1;i++) { 
			postMC(DSR[i],DY[i+1],solver);
			postMCInv(DY[i+1],DSR[i],solver);
		}
		for (int i=0; i<n-1; i++) {
			postSR(DeltaSY[i],DSR[i], solver);
		}

		for(int i=0;i<n;i++) {
			postARK(DY[i],DK[i], DX[i], solver);
		}

		// int cpt=0;
		// for (int J=0;J<n*BC;J++) {
		// 	r=J/BC;
		// 	j=J%BC;
		// 	indP=(J/KC);
		// 	for (int i=0;i<4;i++) {		
		// 		if (J<BC*n-1) {
		// 			if (J%KC==(KC-1)) {
		// 				LCF.ifThenElse(
		// 						ICF.arithm(DK[r][(i+1)%4][j],"=",0),
		// 						LCF.and(
		// 								ICF.arithm(pk[indP][i], "=", 0),
		// 								ICF.arithm(DeltaSK[indP][i],"=",0)
		// 								),
		// 								ICF.table(new IntVar[]{DK[r][(i+1)%4][j], DeltaSK[indP][i],pk[indP][i]},tupleSB, strategy)
		// 						);
		// 			}
		// 		}
		// 		if (J%KC==0 && J/KC>0) {
		// 			postXorByte(DK[(J-KC)/BC][i][(J-KC)%BC],DeltaSK[J/KC-1][i],DK[r][i][j],solver);
		// 		}

		// 		if (J/KC>0 && J%KC>0) {
		// 			postXorByte(DK[(J-1)/BC][i][(J-1+BC)%BC],DK[(J-KC)/BC][i][(J-KC+BC)%BC],DK[r][i][j],solver);
		// 		}
		// 	}	
		// }

		for(r = 0; r < n-1; r++){
			for(int i = 0; i < 4; i++){
				for(j = 0; j < 4; j++){
					int pij = Pk[4*j+i];
					int j2 = pij/4;
					int i2 = pij%4;
					solver.post(IntConstraintFactory.arithm(DK[r+1][i2][j2], "=", DK[r][i][j]));
				}
			}
		}

		for (r=0;r<n;r++) {
			for (int i=0;i<4;i++) {
				for (j=0;j<BC;j++) {
					allProbas[cproba++]=p[r][i][j];
					// if ((BC*(r)+j)%KC==KC-1) {
					// 	allProbas[cproba++]=pk[(BC*(r)+j)/KC][i];
					// }
				}	
			}
		}
		obj=VF.bounded("objs2", 0, 2*nbSB, solver);
		solver.post(ICF.sum(allProbas, "=",obj));
		solver.post(IntConstraintFactory.arithm(obj, ">=", 8*nbSB-128)); 
	}




	@Override
	public void configureSearch() {

		solver.set(IntStrategyFactory.domOverWDeg(solver.retrieveIntVars(), System.currentTimeMillis()));//,IntStrategyFactory.minDom_UB(getAllVar()));
		solver.set(IntStrategyFactory.lastKConflicts(solver,12,solver.getStrategy()));	
		this.level=Level.SILENT;

	}



	@Override
	public void createSolver() {
		solver=new Solver("step2");
	}



	@Override
	public void prettyOut() {
		for(Solution s:solver.getSolutionRecorder().getSolutions()){
			int pr=s.getIntVal(obj)-(8*nbSB);
			System.out.println("Proba: "+pr);

			writer.println("Proba :2^"+pr);
			writer.println("Round 0:");
			writer.println("       deltaX[0]                deltaK[0]                deltaY[0]                DeltaSY[0]                pSY[0] ");
			for (int j=0; j<4; j++){
				for (int k=0; k<BC; k++) writer.printf("%5s",s.getIntVal(DY[0][j][k])+" "); writer.print("     ");
				for (int k=0; k<BC; k++) writer.printf("%5s",s.getIntVal(DK[0][j][k])+" "); 

				writer.print("     ");
				//writer.printf("%5s",s.getIntVal(pk[0][j])+" "); writer.print("     ");
				//writer.print("sk="+DeltaSK[0][j].getValue()+" ");
				for (int k=0; k<BC; k++) writer.printf("%5s",s.getIntVal(DX[0][j][k])+" "); writer.print("     ");
				for (int k=0; k<BC; k++) writer.printf("%5s",s.getIntVal(DeltaSY[0][j][k])+" "); writer.print("     ");
				for (int k=0; k<BC; k++) writer.printf("%5s",s.getIntVal(p[0][j][k])+" "); writer.println("     ");

			}
			for (int i=1; i<n; i++){
				writer.println("Round "+i+":");
				writer.println("       deltaSR["+(i-1)+"]                deltaX["+(i)+"]                deltaK["+i+"]                deltaY["+i+"]                deltaSY["+i+"]                 pSY["+i+"]");

				for (int j=0; j<4; j++){
					for (int k=0; k<BC; k++) writer.printf("%5s",s.getIntVal(DSR[i-1][j][k])+" "); writer.print("      ");
					for (int k=0; k<BC; k++) writer.printf("%5s",s.getIntVal(DY[i][j][k])+" "); writer.print("      ");
					for (int k=0; k<BC; k++) writer.printf("%5s",s.getIntVal(DK[i][j][k])+" "); writer.print("      ");

					for (int k=0; k<BC; k++) writer.printf("%5s",s.getIntVal(DX[i][j][k])+" "); 
					writer.print("      ");
					if (i<n-1)
						for (int k=0; k<BC; k++) writer.printf("%5s",s.getIntVal(DeltaSY[i][j][k])+" "); writer.print("      ");
						for (int k=0; k<BC; k++) writer.printf("%5s",s.getIntVal(p[i][j][k])+" "); 

						writer.println();
				}
			}

			writer.flush();	
		}
	}



	@Override
	public void solve() {
		System.out.println("Beginning solving");
		boolean b = solver.findSolution();
		if(!b){ 
			if(solver.hasReachedLimit()){
				System.out.println("No solution found, but a limit has been reached");
			}
			System.out.println("No solutions !!");
			Chatterbox.printShortStatistics(solver);
		}
		else{
			prettyOut();
			Chatterbox.printShortStatistics(solver);
		}
		//solver.findOptimalSolution(ResolutionPolicy.MAXIMIZE, obj);
//		solver.plugMonitor(new IMonitorSolution() {
//            @Override
//            public void onSolution() {
//            	prettyOut();
//            	//System.out.println("solution");
//            }
//        });
//		boolean b=solver.findSolution();
//		int cpt=0;
//		
//		while (b)
//			b=solver.nextSolution();
		//System.out.println();
		//prettyOut();
		//System.out.println("\n Total time posting xor byte : "+ 		totalPXB+ "s");
		//Chatterbox.printShortStatistics(solver);


	}




}

