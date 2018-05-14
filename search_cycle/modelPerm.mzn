include "value_precede_chain.mzn";
array[0..15] of var 0..15: Pk; %Permutation (key schedule)
array[0..15] of var 0..15: Ps; %Permutation (instead of ShiftRows)


int: n; % Number of rounds
%int: objective = 14; % Max objective value
int: objective_bound;
int: KEY_BITS=128; % Number of bits in the key

int: BLOCK_BITS=128; % Number of bits in the blocks
int: KC=KEY_BITS div 32; % Number of columns per round of key schedule
int: BC=BLOCK_BITS div 32; % Number of columns per round
int: NBK =  KC + n-1;  % Maximum number of components per byte of the key
int: NK =  n; % Number of columns with SBoxes in the key

array[0..n-1, 0..15] of var 0..15: Ks; %Id, Pk, Pk^2, ...
constraint forall (i in 0..15) (Ks[0,i] = i);
constraint forall (r in 1..n-1, i in 0..15) (Ks[r,Pk[i]] = Ks[r-1,i]);


array[0..n-2, 0..BC-1, 0..3] of var 0..1: deltaY; % State before ARK
array[0..n-1, 0..BC-1,0..3] of var 0..1: deltaX; % State after ARK
array[0..n-1, 0..BC-1,0..3] of var 0..1: deltaSR; % State after ShiftRows

array[0..BC-1,0..3] of var 0..1: deltaK; % Key
array[0..127] of var 0..1: deltaK_bits;

constraint forall (i in 0..3, j in 0..3) (
  (deltaK[j,i] == 0) <-> ((deltaK_bits[8*(4*j+i)] == 0) /\ (deltaK_bits[8*(4*j+i)+1] == 0) /\ (deltaK_bits[8*(4*j+i)+2] == 0)
            /\ (deltaK_bits[8*(4*j+i)+3] == 0) /\ (deltaK_bits[8*(4*j+i)+4] == 0) /\ (deltaK_bits[8*(4*j+i)+5] == 0)
            /\ (deltaK_bits[8*(4*j+i)+6] == 0) /\ (deltaK_bits[8*(4*j+i)+7] == 0))
);


array[0..n-1,0..BC-1] of var 0..4: colX; %   colX[r][j] = The sum for i in 0..3 of DeltaX[r][i][j]
array[0..n-1,0..BC-1] of var 0..4: colSRX; % colSRX[r][j] = The sum for i in 0..3 of SR(DeltaX)[r][i][j]
array[0..n-2,0..BC-1] of var 0..4: eqs; % colSRX[r][j] = The sum for i in 0..3 of SR(DeltaX)[r][i][j]

%The objective variable
var 1..KEY_BITS div 6 +1+1 :obj;
var 1..16: wK; %number of actives bytes in the key

%The sum to minimize (non zero values going inside Sboxes)
constraint obj = sum(r in 0..n-1, j in 0..BC-1) (colSRX[r,j]);
%constraint obj=objective;
constraint obj <= 19;
constraint wK = sum(i in 0..3, j in 0..3) (deltaK[j,i]);

% regular XOR, except 1 + 1 can be either 0 or 1.
predicate XOR(var 0..1: a,var 0..1: b,var 0..1: c) =
      a+b+c!=1
    ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation of the redundant variables.
constraint forall (r in 0..n-1, j in 0..3) (
                    colX[r,j]=sum(i in 0..3)(deltaX[r,j,i]) /\
                    colSRX[r,j]=sum(i in 0..3)(deltaSR[r,j,i])
);

%%%%%%%%%%%%%%%%ARK
constraint forall(r in 1..n-1,j in 0..3, i in 0..3) (
    let {var 0..15 : pij = Ks[r,4*j+i], var 0..3: j1 = pij div 4, var 0..3: i1 = pij mod 4} in
       XOR(deltaY[r-1,j,i], deltaK[j1,i1], deltaX[r,j,i])
)
;

%%%%%%%%%%%%%%% Ps (instead of ShiftRows)
% constraint  forall(r in 0..n-1,j in 0..3, i in 0..3) (
%        deltaSR[r,j,i]=deltaX[r,((j+i) mod BC),i]
% );
constraint forall (r in 1..n-1, i in 0..3, j in 0..3) (
  let {var 0..15: pij = Ps[4*j+i], var 0..3: j1 = pij div 4, var 0..3: i1 = pij mod 4} in
    deltaSR[r,j1,i1] = deltaX[r,j,i]
);

%%%%%%%%%%%%%%% Mixcolumns
constraint forall(r in 0..n-2, j in 0..3) (
          (colSRX[r,j] + sum(i in 0..3)(deltaY[r,j,i])) in {0,5,6,7,8}
  )
;



constraint forall(r in 0..n-2, j in 0..3) (

(((deltaSR[r,j,0] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+2]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+2]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+3]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+3]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+3]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+3]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+3]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+4] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0)))
 /\
  (((deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j+1]] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j+1]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j+1]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j+1]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+1]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j+2]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+2]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j+1]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j+3]] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j+3]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j+3]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j+3]+3] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+3]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+2]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0)))
 /\
  (((deltaSR[r,j,0] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+3]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j+1]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+1]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+1]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+1]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+2]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,0] == 1) /\ (deltaX[r+1,j,0] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]] + deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+4] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+1] + deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+5] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+5] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+2] + deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+6] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j]+3] + deltaK_bits[8*Ks[r,4*j]+4] + deltaK_bits[8*Ks[r,4*j]+7] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+2]+1] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+2]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+2]+5] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+2]+6] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,1] == 1) /\ (deltaX[r+1,j,1] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]] + deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+1] + deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+2] + deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+6] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+3] + deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+1]+4] + deltaK_bits[8*Ks[r,4*j+1]+5] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0)))
 /\
  (((deltaSR[r,j,1] == 0) /\ (deltaSR[r,j,2] == 0) /\ (deltaSR[r,j,3] == 0) /\ (deltaY[r,j,2] == 1) /\ (deltaX[r+1,j,2] == 0) /\ (deltaY[r,j,3] == 1) /\ (deltaX[r+1,j,3] == 0)) -> (
(deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]] + deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+1] + deltaK_bits[8*Ks[r,4*j+3]+6] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+1] + deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+2] + deltaK_bits[8*Ks[r,4*j+3]+6] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+2] + deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+6] + deltaK_bits[8*Ks[r,4*j+3]+3] + deltaK_bits[8*Ks[r,4*j+3]+7] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+3] + deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+7] + deltaK_bits[8*Ks[r,4*j+3]+4] mod 2 == 0) /\
(deltaK_bits[8*Ks[r,4*j+2]+4] + deltaK_bits[8*Ks[r,4*j+2]+5] + deltaK_bits[8*Ks[r,4*j+3]+5] mod 2 == 0)))

);


% The search heuristic. We first assign the variables corresponding to the SBoxes positions, and then the rest, always starting by 0 (indomain_min)

solve :: seq_search([
int_search(

          [colSRX[r,j] | r in 0..n-1, j in 0..BC-1] , smallest, indomain_min, complete) ,
         int_search(
          array1d(deltaSR)++
          [deltaK[J mod BC,i] | J in 0..BC*n -1, i in 0..3 where J mod KC == KC-1] ++
          [colX[J div BC, J mod BC] | J in 0..BC*n -1 where J mod KC==KC-1] ++
          [deltaY[J div BC ,J mod BC,i] | J in 0..BC*(n-1) -1, i in 0..3 where J mod KC == KC-1] ++
          [deltaK[J mod BC,i] | J in 0..BC*n -1, i in 0..3 where J mod KC != KC-1] ++
          [colX[J div BC, J mod BC] | J in 0..BC*n -1 where J mod KC!=KC-1] ++
          [deltaX[r,j,i] | r in 0..n-1,  j in 0..BC-1, i in 0..3]  ++
          [deltaY[J div BC ,J mod BC,i ] | J in 0..BC*(n-1) -1, i in 0..3 where J mod KC != KC-1] , smallest
, indomain_min, complete)])

satisfy;
%minimize obj;

%output ["deltaY=array3d(0.."++show(n-2)++", 0..3, 0..3,"++show(deltaY)++");\ndeltaX=array3d(0.."++show(n-1)++", 0..3, 0..3,"++show(deltaX)++");\ndeltaK=array3d(0.."++show(n-1)++", 0..3, 0..3,"++show(deltaK)++");"];

%Full output
%output [show([if r>0 then deltaY[r-1,j,i] else 0 endif | j in 0..BC-1] )  ++ "   "++show([deltaK[Ks[r,4*j+i] div 4, Ks[r,4*j+i] mod 4] | j in 0..BC-1]) ++ "   "++ show([deltaX[r,j,i]| j in 0..BC-1]) ++ "   " ++ show([deltaSR[r,j,i] | j in 0..BC-1] ) ++ "   " ++  if i mod 3==0 /\ i>0 then "\n\n" else "\n" endif  | r in 0..n-1, i in 0..3] ++ [show(obj) ++ "\nkey : " ++ show(wK) ++ "\n"];

%output [show(obj)++"\n"];
output [show(wK)++"\n"];
%output [ifshow(deltaK[0,j,i]) | j in 0..3 ,i in 0..3] ++ ["\n"];
%output [if 1 < n then show(deltaK[1,j,i]) else "" endif | j in 0..3 ,i in 0..3] ++ ["\n"];
%output [if 2 < n then show(deltaK[2,j,i]) else "" endif | j in 0..3 ,i in 0..3] ++ ["\n"];
%output [if 3 < n then show(deltaK[3,j,i]) else "" endif | j in 0..3 ,i in 0..3] ++ ["\n"];
%output [if 4 < n then show(deltaK[4,j,i]) else "" endif | j in 0..3 ,i in 0..3] ++ ["\n"];
%output [if 5 < n then show(deltaK[5,j,i]) else "" endif | j in 0..3 ,i in 0..3] ++ ["\n"];
%output [if 6 < n then show(deltaK[6,j,i]) else "" endif | j in 0..3 ,i in 0..3] ++ ["\n"];

%output [if 0 < n then show(deltaSR[0,j,i]) else "" endif | j in 0..3 ,i in 0..3] ++ ["\n"];
%output [if 1 < n then show(deltaSR[1,j,i]) else "" endif | j in 0..3 ,i in 0..3] ++ ["\n"];
%output [if 2 < n then show(deltaSR[2,j,i]) else "" endif | j in 0..3 ,i in 0..3] ++ ["\n"];
%output [if 3 < n then show(deltaSR[3,j,i]) else "" endif | j in 0..3 ,i in 0..3] ++ ["\n"];
%output [if 4 < n then show(deltaSR[4,j,i]) else "" endif | j in 0..3 ,i in 0..3] ++ ["\n"];
%output [if 5 < n then show(deltaSR[5,j,i]) else "" endif | j in 0..3 ,i in 0..3] ++ ["\n"];
%output [if 6 < n then show(deltaSR[6,j,i]) else "" endif | j in 0..3 ,i in 0..3] ++ ["\n"];
%output [if 7 < n then show(deltaSR[7,j,i]) else "" endif | j in 0..3 ,i in 0..3] ++ ["\n"];
