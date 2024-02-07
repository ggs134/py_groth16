#q
curve_order = 21888242871839275222246405745257275088548364400416034343698204186575808495617
q = curve_order
#p
p = 21888242871839275222246405745257275088696311157297823662689037894645226208583


Ap = [
    [-60.0, 110.0, -60.0, 10.0],
    [96.0, -136.0, 60.0, -8.0],
    [0.0, 0.0, 0.0, 0.0],
    [-72.0, 114.0, -48.0, 6.0],
    [48.0, -84.0, 42.0, -6.0],
    [-12.0, 22.0, -12.0, 2.0]
]

Bp = [
    [36.0, -62.0, 30.0, -4.0],
    [-24.0, 62.0, -30.0, 4.0],
    [0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0]
]

Cp = [
    [0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0],
    [-144.0, 264.0, -144.0, 24.0],
    [576.0, -624.0, 216.0, -24.0],
    [-864.0, 1368.0, -576.0, 72.0],
    [576.0, -1008.0, 504.0, -72.0]
]

Z = [3456.0, -7200.0, 5040.0, -1440.0, 144.0]
R = [1, 3, 35, 9, 27, 30]

Ap = [ [int(num) % curve_order for num in vec] for vec in Ap ]
Bp = [ [int(num) % curve_order for num in vec] for vec in Bp ]
Cp = [ [int(num) % curve_order for num in vec] for vec in Cp ]
Zp = [ int(num) % curve_order for num in Z ]
Rp = [ int(num) % curve_order for num in R ]

Z = IntegerModRing(q) #Z mod q
R.<x> = PowerSeriesRing(Z) #[a1, a2, a3 ...] = a0x^0 + a1x^2 + a2x^3 ...

Ax = matrix(Z, Ap)
Bx = matrix(Z, Bp)
Cx = matrix(Z, Cp)
Zx = vector(Z, Zp)
Rx = vector(Z, Rp)

Rax = R(list(Rx*Ax))
Rbx = R(list(Rx*Bx))
Rcx = R(list(Rx*Cx))
Px = Rax * Rbx - Rcx

px = Px.polynomial()
zx = R(Zp).polynomial()

Hx = px // zx # 12*x^3 + 4836*x^2 + 384*x + 720
# print(Hx.list())
remainder = px % zx # 0

print("The remainder should be 0 : {}".format(remainder == 0))

############################
### 1.2 CRS CONSTRUCTION ###
############################

alpha = Z(3926)
beta = Z(3604)
gamma = Z(2971)
delta = Z(1357)
x_val = Z(3721)

# alpha = random.randint(1,q-2)
# beta = random.randint(1,q-2)
# gamma = random.randint(1,q-2)
# delta = random.randint(1,q-2)
# x_val = random.randint(1,q-2)

tau = [alpha, beta, gamma, delta, x_val]

# tau : [3351, 1339, 4343, 2771, 380]
# tau : [3806, 1193, 1247, 1027, 1707]
# tau : [1203, 320, 4771, 4621, 183]

print("tau : {}".format(tau))

Ax_val = []
Bx_val = []
Cx_val = []
# Zx_val
# Hx_val

for ax in Ax.rows():
    ax_single = R(list(ax))(x_val)
    Ax_val.append(ax_single)

for bx in Bx.rows():
    bx_single = R(list(bx))(x_val)
    Bx_val.append(bx_single)

for cx in Cx.rows():
    cx_single = R(list(cx))(x_val)
    Cx_val.append(cx_single)

Zx_val = R(Zx.list())(x_val)
Hx_val = Hx(x_val)

numGates = len(Ax.columns())
numWires = len(Ax.rows())

print("numGates : {}".format(numGates))
print("numWires : {}".format(numWires))

# numGates : 5
# numWires : 6

# print("Ax_val : {}".format(Ax_val))
# print(Bx_val)
# print(Cx_val)
# print(Zx_val)
# print(Hx_val)

# Elliptic curve definition G1, G2
# % k=2;
# % % Elliptic curve: y^2=x^3+a*x+0
# % a=1;
# % b=0;
# % p=71; % Field size for elliptic curve arithmetic. Choose any prime p such that mod(p,4)==3.

# This is BN128 curve
# field_modulus = 21888242871839275222246405745257275088696311157297823662689037894645226208583
# y² = x³ + 3 (mod field_modulus)
EC = EllipticCurve(GF(p), [0,0,0,0,3])
g = EC(1,2)

# g = EC.points()[12] #(11 : 8 : 1), order is 72, which is max among points

F2 = GF(p^2,"i",modulus=x^2 + 1)
TwistB = 3*F2("9+i")^(-1)
G2 = EllipticCurve(F2,[0,TwistB])

P2x = F2("11559732032986387107991004021392285783925812861821192530917403151452391805634*i + 10857046999023057135944570762232829481370756359578518086990519993285655852781")
P2y = F2("4082367875863433681332203403145435568316851327593401208105741076214120093531*i + 8495653923123431417604973247489272438418190587263600148770280649306958101930")
g2 = G2(P2x,P2y)
h = g2

# F.<z> = GF(p^2, modulus = x^2 + 1)
# ECExt = EllipticCurve(F,[1,0])
# ECExt = EC.base_extend(F)
# ECExt = EC.change_ring(F)
# maxpoints = [p for p in ECExt.points() if p.order() == 72]
# h = [mp for mp in maxpoints if mp[0] == 60][0] #extended generator

# ECExt2 = EllipticCurve(F,[1,0], gen=h)

pointInf1 = g * curve_order
pointInf2 = h * curve_order

### sigma1_1 ###
sigma1_1 = [g*alpha, g*beta, g*delta]

sigma1_2 = []
sigma1_3 = []
sigma1_4 = []
sigma1_5 = []

sigma2_1 = [beta*h, gamma*h, delta*h]
sigma2_2 = []

### sigma1_2 ###
for i in range(numGates):
    val = x_val^i
    sigma1_2.append(val * g)

# print(sigma1_2)

### sigma1_3 ###
VAL = [0]*numWires
for i in range(numWires):
    if i in [0, numWires-1]:
        val = (beta*Ax_val[i] + alpha*Bx_val[i] + Cx_val[i]) // gamma
        VAL[i] = val
        sigma1_3.append(val * g)
    else:
        sigma1_3.append(0)

# print("sigma1_3 : {}".format(sigma1_3))
# print(sigma1_3)

#### sigma1_4 ###
for i in range(numWires):
    if i in [0, numWires-1]:
        sigma1_4.append(0)
    else:
        val = (beta*Ax_val[i] + alpha*Bx_val[i] + Cx_val[i]) // delta
        sigma1_4.append(val * g)


# print("sigma1_4 : {}".format(sigma1_4))

#### sigma1_5 ###
for i in range(numGates-1):
    sigma1_5.append(g*(x_val^i * Zx_val // delta))

# print(sigma1_5)

#sigma2-2
for i in range(numGates):
    # sigma2_2.append(h*(Z(x_val^i)))
    sigma2_2.append(h*(x_val^i))

print("CRS(proving / verification key) : ")
print("Sigma1_1 : {}".format(sigma1_1))
print("Sigma1_2 : {}".format(sigma1_2))
print("Sigma1_3 : {}".format(sigma1_3))
print("Sigma1_4 : {}".format(sigma1_4))
print("Sigma1_5 : {}".format(sigma1_5))
print("Sigma2_1 : {}".format(sigma2_1))
print("Sigma2_2 : {}".format(sigma2_2))

# Sigma1_1 : [(24 : 28 : 1), (10 : 67 : 1), (39 : 12 : 1)]
# Sigma1_2 : [(11 : 8 : 1), (23 : 64 : 1), (51 : 28 : 1), (11 : 8 : 1), (23 : 64 : 1)]
# Sigma1_3 : [(0 : 1 : 0), 0, 0, 0, 0, (0 : 1 : 0)]
# Sigma1_4 : [0, (0 : 1 : 0), (0 : 1 : 0), (0 : 1 : 0), (0 : 1 : 0), 0]
# Sigma1_5 : [(0 : 1 : 0), (0 : 1 : 0), (0 : 1 : 0), (0 : 1 : 0)]
# Sigma2_1 : [(61 : 67*z : 1), (9 : 55*z : 1), (32 : 12*z : 1)]
# Sigma2_2 : [(60 : 8*z : 1), (48 : 64*z : 1), (20 : 28*z : 1), (60 : 8*z : 1)]

##############################
### 1.3 CRS VALIDITY CHECK ###
##############################

# Ax_val_vec = vector(ZZ, Ax_val)
# Bx_val_vec = vector(ZZ, Bx_val)
# Cx_val_vec = vector(ZZ, Cx_val)

# Zx_val = Z(Zx_val)
# Hx_val = Z(Hx_val)

Ax_val = vector(Z, Ax_val)
Bx_val = vector(Z, Bx_val)
Cx_val = vector(Z, Cx_val)

Zx_val = Z(Zx_val)
Hx_val = Z(Hx_val)

# print("Zx_val : {}".format(Zx_val))
# print("Hx_val : {}".format(Hx_val))

# print("Rx*Ax_val_vec : {}".format(Rx*Ax_val_vec))
# print("Rx*Bx_val_vec : {}".format(Rx*Bx_val_vec))
# print("Cx*Cx_val_vec : {}".format(Rx*Cx_val_vec))

# lhs = Z((Rx*Ax_val_vec)*(Rx*Bx_val_vec)-(Rx*Cx_val_vec))
# rhs = Zx_val*Hx_val

lhs = Z((Rx*Ax_val)*(Rx*Bx_val)-(Rx*Cx_val))
rhs = Zx_val*Hx_val

# print("lhs : {}".format(lhs))
# print("rhs : {}".format(rhs))

# print("lhs : {}".format(lhs))
# print("rhs : {}".format(rhs))

if lhs == rhs:
    print('CRS is valid (Setup successful)')
    print()
else:
    print('CRS is invalid (Setup fails)')
    print()

##################
### 2. PROVING ###
##################

#random number created by proover
r = Z(4106)
s = Z(4565)

# r = random.randint(0, q-1)
# s = random.randint(0, q-1)

print("r, s : {}, {}".format(r, s))

#Build Proof_A
proof_A = sigma1_1[0]

for i in range(numWires):
    temp = pointInf1    
    for j in range(numGates):
        temp = temp + (Ax[i][j] * sigma1_2[j])    
    proof_A = proof_A + (Rx[i] * temp)
proof_A = proof_A + (r * sigma1_1[2])

#Build proof_B
proof_B = sigma2_1[0]
for i in range(numWires):
    temp = pointInf2   
    for j in range(numGates):
        temp = temp + (Bx[i][j] * sigma2_2[j])    
    proof_B = proof_B + (Rx[i] * temp)
proof_B = proof_B + (s * sigma2_1[2])

#Build temp_proof_B
temp_proof_B = sigma1_1[1]
for i in range(numWires):
    temp = pointInf1
    for j in range(numGates):
        temp = temp + (Bx[i][j] * sigma1_2[j])    
    temp_proof_B = temp_proof_B + (Rx[i] * temp)
temp_proof_B = temp_proof_B + (s * sigma1_1[2])

#Build proof_C
proof_C = (s * proof_A) + (r * temp_proof_B) - (r*s*sigma1_1[2])

for i in range(1,numGates-1):
    proof_C = proof_C + (Rx[i] * sigma1_4[i])

for i in range(numGates-1):
    proof_C = proof_C + (Hx[i] * sigma1_5[i])

proof = [proof_A, proof_B, proof_C]

print("Proofs [proof_A, proof_B, proof_C] : ")
print(proof)
print()

####################################
### 2.1 PROOF COMPLETENESS CHECK ###
####################################

A = alpha + Rx*Ax_val + r*delta
B = beta + Rx*Bx_val + s*delta

# C = 1/delta * Rx[1:numWires-1]*(beta*Ax_val[1:numWires-1] + alpha*Bx_val[1:numWires-1] + Cx_val[1:numWires-1]) + Hx_val*Zx_val + A*s + B*r + (-r*s*delta)
C = 1/delta * \
    (Rx[1:numWires-1] * \
        (beta*Ax_val[1:numWires-1] + alpha*Bx_val[1:numWires-1] + Cx_val[1:numWires-1]) \
        + Hx_val*Zx_val \
    ) \
     + A*s \
     + B*r \
     + (-r*s*delta)

lhs = A*B

rhs = alpha*beta #14149304

rpub = [Rx[0], Rx[-1]]
valpub = [VAL[0], VAL[-1]]

rhs = rhs + gamma*vector(rpub)*vector(valpub)  #12058091336480024
rhs = rhs + C*delta #21888242871839275222246405745257275088548364400416033032405666501928354297837

result = (proof_A == g*A) and (proof_B == B*h) and (proof_C == C*g)

print("proof completeness check : {}".format(result and lhs==rhs))

#TODO : not working

######################
##### 3. VERIFY ######
######################

# def weil(point1, point2):
#     val = EC.isomorphism_to(G2)(point1).weil_pairing(point2, p)
#     return val

# LHS = weil(proof_A, proof_B)
# RHS = weil(sigma1_1[0], sigma2_1[0])

# temp = pointInf1

# for i in [0, numWires-1]:
#   temp = temp + Rx[i]*sigma1_3[i]

# RHS = RHS * weil(temp, sigma2_1[1])
# RHS = RHS * weil(proof[2], sigma2_1[2])

# print(LHS)
# print(RHS)

# print("Verification result (RHS == LHS)? : {}".format(RHS == LHS))