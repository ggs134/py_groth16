from py_ecc.bn128 import multiply, G1, G2, add, pairing, neg, curve_order, final_exponentiate
import galois
import numpy as np

#Taking some time..
GF = galois.GF(curve_order)

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

Ax = [ [int(num) % curve_order for num in vec] for vec in Ap ]
Bx = [ [int(num) % curve_order for num in vec] for vec in Bp ]
Cx = [ [int(num) % curve_order for num in vec] for vec in Cp ]
Zx = [ int(num) % curve_order for num in Z ]
Rx = [ int(num) % curve_order for num in R ]

npAx = GF(np.array(Ax))
npBx = GF(np.array(Bx))
npCx = GF(np.array(Cx))
npZx = GF(np.array(Zx))
npRx = GF(np.array(Rx))

npRax = npAx.transpose().dot(npRx)
npRbx = npBx.transpose().dot(npRx)
npRcx = npCx.transpose().dot(npRx)

Rax = galois.Poly(npRax, order="asc")
Rbx = galois.Poly(npRbx, order="asc")
Rcx = galois.Poly(npRcx, order="asc")

Zx = galois.Poly(npZx, order="asc")
Rx = galois.Poly(npRx, order="asc")

Px = Rax * Rbx - Rcx

Hx = Px // Zx       #quotient
Remainder = Px % Zx #remainder

print("Px % Zx  = 0 ?  {}".format(Remainder == 0))

############################
### 1.2 CRS CONSTRUCTION ###
############################

alpha = GF(3926)
beta = GF(3604)
gamma = GF(2971)
delta = GF(1357)
x_val = GF(3721)

tau = [alpha, beta, gamma, delta, x_val]

#### DONE SO FAR ####

Ax_val = []
Bx_val = []
Cx_val = []

for i in range(len(Ax)):
    ax_single = galois.Poly(Ax[i], order="asc", field=GF)(x_val)
    Ax_val.append(ax_single)

for i in range(len(Bx)):
    bx_single = galois.Poly(Bx[i], order="asc", field=GF)(x_val)
    Bx_val.append(bx_single)

for i in range(len(Cx)):
    cx_single = galois.Poly(Cx[i], order="asc", field=GF)(x_val)
    Cx_val.append(cx_single)

Zx_val = Zx(x_val)
Hx_val = Hx(x_val)

numGates = len(Ax[0])
numWires = len(Ax)

sigma1_1 = [multiply(G1, int(alpha)), multiply(G1, int(beta)), multiply(G1, int(delta))]
sigma1_2 = []
sigma1_3 = []
sigma1_4 = []
sigma1_5 = []

sigma2_1 = [multiply(G2, int(alpha)), multiply(G2, int(beta)), multiply(G2, int(delta))]
sigma2_2 = []

#sigma1_2
for i in range(numGates):
    val = x_val ** i
    sigma1_2.append(multiply(G1, int(val)))

#sigma1_3
VAL = [GF(0)]*numWires
for i in range(numWires):
    if i in [0, numWires-1]:
        val = (beta*Ax_val[i] + alpha*Bx_val[i] + Cx_val[i]) / gamma
        VAL[i] = val
        sigma1_3.append(multiply(G1, int(val)))
    else:
        sigma1_3.append((0, 0))

#sigma1_4
for i in range(numWires):
    if i in [0, numWires-1]:
        sigma1_4.append((0, 0))
    else:
        val = (beta*Ax_val[i] + alpha*Bx_val[i] + Cx_val[i]) / delta
        sigma1_4.append(multiply(G1, int(val)))

#sigma1_5
for i in range(numGates-1):
    sigma1_5.append(multiply(G1, int((x_val**i * Zx_val) / delta)))

#sigma2-2
for i in range(numGates):
    # sigma2_2.append(h*(Z(x_val^i)))
    sigma2_2.append(multiply(G2, int(x_val**i)))

lhs = Rax(x_val)*Rbx(x_val) - Rcx(x_val)
rhs = Zx_val*Hx_val

print("Is lhs == rhs ? : {}".format(lhs == rhs))

### 2. PROVING ###
# TODO : imp

### 2.1 PROOF COMPLETENESS CHECK ###
# TODO : imp

##### 3. VERIFY ######
# TODO : imp




