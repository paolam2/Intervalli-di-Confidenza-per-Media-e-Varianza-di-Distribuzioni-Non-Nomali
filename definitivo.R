# Simulazione per Tesi di Laurea: 
# INTERVALLI DI CONFIDENZA PER MEDIA E VARIANZA DI DISTRIBUZIONI NON NORMALI

# Autore adattamento: 
# PAOLA MARIA LEPORE

# Basato su codice originale di: 
# Prof. JOSE DIAS CURTO, Instituto Universitario de Lisboa (ISCTE-IUL), 
# BRU-UNIDE, Lisboa, Portugal - Università di [Nome]

# Articolo del prof. Curto disponibile al link:
# https://doi.org/10.1080/03610918.2021.1963448

################################################################################
# VARIANZA
# TABELLA 1: DISTORSIONE DEGLI STIMATORI
M = 10000
n = 100

# NORMALE(0,1), CURTOSI = 3
sim1 = matrix(NA, M)
sim2 = matrix(NA, M)
sim3 = matrix(NA, M)
for(i in 1:M){
  camp = rnorm(n, 0, 1)
  perc = 1/(2*(n-4)^(1/2))
  gamma41 = n*sum((camp-mean(camp))^4)/(sum((camp-mean(camp))^2))^2
  gamma42 = n*sum((camp-mean(camp, perc, na.rm = FALSE))^4)/(sum((camp-mean(camp))^2))^2
  gamma43 = n*sum((camp-median(camp))^4)/(sum((camp-mean(camp))^2))^2
  sim1[i] = gamma41 - 3
  sim2[i] = gamma42 - 3
  sim3[i] = gamma43 - 3
}
out1 = mean(sim1)
out2 = mean(sim2)
out3 = mean(sim3)

# UNIFORME(0,1), CURTOSI = 9/5 = 1.8
sim1u = matrix(NA, M)
sim2u = matrix(NA, M)
sim3u = matrix(NA, M)
for(i in 1:M){
  camp = runif(n)
  perc = 1/(2*(n-4)^(1/2))
  gamma41u = n*sum((camp-mean(camp))^4)/(sum((camp-mean(camp))^2))^2
  gamma42u = n*sum((camp-mean(camp, perc, na.rm = FALSE))^4)/(sum((camp-mean(camp))^2))^2
  gamma43u = n*sum((camp-median(camp))^4)/(sum((camp-mean(camp))^2))^2
  sim1u[i] = gamma41u - 1.8
  sim2u[i] = gamma42u - 1.8
  sim3u[i] = gamma43u - 1.8
}
out1u = mean(sim1u)
out2u = mean(sim2u)
out3u = mean(sim3u)

# CHI-QUADRO(1), CURTOSI = 3+12/1 = 15
sim1c = matrix(NA, M)
sim2c = matrix(NA, M)
sim3c = matrix(NA, M)
for(i in 1:M){
  camp = rchisq(n,1)
  perc = 1/(2*(n-4)^(1/2))
  gamma41c = n*sum((camp-mean(camp))^4)/(sum((camp-mean(camp))^2))^2
  gamma42c = n*sum((camp-mean(camp, perc, na.rm = FALSE))^4)/(sum((camp-mean(camp))^2))^2
  gamma43c = n*sum((camp-median(camp))^4)/(sum((camp-mean(camp))^2))^2
  sim1c[i] = gamma41c - 15
  sim2c[i] = gamma42c - 15
  sim3c[i] = gamma43c - 15
}
out1c = mean(sim1c)
out2c = mean(sim2c)
out3c = mean(sim3c)

#  T-STUDENT(5), CURTOSI = (3*(v-2))/(v-4)
v = 5
sim1t = matrix(NA, M)
sim2t = matrix(NA, M)
sim3t = matrix(NA, M)
for(i in 1:M){
  camp = rt(n,v)
  perc = 1/(2*(n-4)^(1/2))
  gamma41t = n*sum((camp-mean(camp))^4)/(sum((camp-mean(camp))^2))^2
  gamma42t = n*sum((camp-mean(camp, perc, na.rm = FALSE))^4)/(sum((camp-mean(camp))^2))^2
  gamma43t = n*sum((camp-median(camp))^4)/(sum((camp-mean(camp))^2))^2
  sim1t[i] = gamma41t - 9
  sim2t[i] = gamma42t - 9
  sim3t[i] = gamma43t - 9
}
out1t = mean(sim1t)
out2t = mean(sim2t)
out3t = mean(sim3t)

# ESPONENZIALE(1), CURTOSI = 9
sim1e = matrix(NA, M)
sim2e = matrix(NA, M)
sim3e = matrix(NA, M)
for(i in 1:M){
  camp = rexp(n)
  perc = 1/(2*(n-4)^(1/2))
  gamma41e = n*sum((camp-mean(camp))^4)/(sum((camp-mean(camp))^2))^2
  gamma42e = n*sum((camp-mean(camp, perc, na.rm = FALSE))^4)/(sum((camp-mean(camp))^2))^2
  gamma43e = n*sum((camp-median(camp))^4)/(sum((camp-mean(camp))^2))^2
  sim1e[i] = gamma41e - 9 
  sim2e[i] = gamma42e - 9
  sim3e[i] = gamma43e - 9
}
out1e = mean(sim1e)
out2e = mean(sim2e)
out3e = mean(sim3e)

# BETA(3,3)
a = 3
b = 3
curtosi = (3*(a+b)*(a+b+1)*(a+1)*(2*b-a))/(a*b*(a+b+2)*(a+b+3))+(a*(a-b)/(a+b))
sim1b = matrix(NA, M)
sim2b = matrix(NA, M)
sim3b = matrix(NA, M)
for(i in 1:M){
  camp = rbeta(n, 3, 3)
  perc = 1/(2*(n-4)^(1/2))
  gamma41 = n*sum((camp-mean(camp))^4)/(sum((camp-mean(camp))^2))^2
  gamma42 = n*sum((camp-mean(camp, perc, na.rm = FALSE))^4)/(sum((camp-mean(camp))^2))^2
  gamma43 = n*sum((camp-median(camp))^4)/(sum((camp-mean(camp))^2))^2
  sim1b[i] = gamma41 - curtosi
  sim2b[i] = gamma42 - curtosi
  sim3b[i] = gamma43 - curtosi
}
out1 = mean(sim1b)
out2 = mean(sim2b)
out3 = mean(sim3b)

################################################################################
# TABELLA 3: STIMA PROBABILITA DI COPERTURA BASATA SUGLI STIMAOTRI DELLA CURTOSI CON SE(1)

n = 10
M = 10000

# N(0,1)
sm0 = matrix(NA, M)
sm1 = matrix(NA, M)
sm2 = matrix(NA, M)
sm3 = matrix(NA, M)
for(i in 1:M){
  camp = rnorm(n, 0, 1)
  var = 1
  perc = 1/(2*(n-4)^(1/2))
  gamma41 = n*sum((camp-mean(camp))^4)/(sum((camp-mean(camp))^2))^2
  gamma42 = n*sum((camp-mean(camp, perc/2, na.rm = F))^4)/(sum((camp-mean(camp))^2))^2
  gamma43 = n*sum((camp-median(camp))^4)/(sum((camp-mean(camp))^2))^2
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  se11 = ((gamma41-(n-3)/n)/(n-1))^(1/2)
  se12 = ((gamma42-(n-3)/n)/(n-1))^(1/2)
  se13 = ((gamma43-(n-3)/n)/(n-1))^(1/2)
  
  LL0 = (n-1)*var(camp)/qchisq((1-alpha/2), n-1)
  UL0 = (n-1)*var(camp)/qchisq((alpha/2), n-1)
  if(LL0 < var & UL0 > var) {
    sm0[i] = 1 } else {
      sm0[i] = 0}
  LL1 = exp(log(var(camp))-z*se11)
  UL1 = exp(log(var(camp)) + z * se11)
  if(LL1 < var & UL1 > var) {
    sm1[i] = 1 } else{
      sm1[i] = 0}
  LL2 = exp(log(var(camp)) - z * se12)
  UL2 = exp(log(var(camp)) + z * se12)
  if(LL2 < var & UL2 > var) {
    sm2[i] = 1} else {
      sm2[i] = 0}
  LL3 = exp(log(var(camp)) - z * se13)
  UL3 = exp(log(var(camp)) + z * se13)
  if(LL3 < var & UL3 > var) {
    sm3[i] = 1} else {
      sm3[i] = 0}
}
out0 = mean(sm0)
out1 = mean(sm1)
out2 = mean(sm2)
out3 = mean(sm3)

# U(0,1)
smu0 = matrix(NA,M)
smu1 = matrix(NA,M)
smu2 = matrix(NA,M)
smu3 = matrix(NA,M)
for(i in 1:M) {
  cu = runif(n, 0, 1)
  var = 1/12
  perc = 1/(2*(n-4)^(1/2))
  gamma41 = n*sum((cu-mean(cu))^4)/(sum((cu-mean(cu))^2))^2
  gamma42 = n*sum((cu-mean(cu, perc, na.rm = F))^4)/(sum((cu-mean(cu))^2))^2
  gamma43 = n*sum((cu-median(cu))^4)/(sum((cu-mean(cu))^2))^2
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  se11 = ((gamma41-(n-3)/n)/(n-1))^(1/2)
  se12 = ((gamma42-(n-3)/n)/(n-1))^(1/2)
  se13 = ((gamma43-(n-3)/n)/(n-1))^(1/2)
  LLu0 = (n-1)*var(cu)/qchisq((1-alpha/2), n-1)
  ULu0 = (n-1)*var(cu)/qchisq((alpha/2), n-1)
  if(LLu0 < var & ULu0 > var) {
    smu0[i] = 1 } else {
      smu0[i] = 0}
  LLu1 = exp(log(var(cu))-z*se11)
  ULu1 = exp(log(var(cu))+z*se11)
  if(LLu1 < var & ULu1 > var) {
    smu1[i] = 1 } else {
      smu1[i] = 0}
  LLu2 = exp(log(var(cu))-z*se12)
  ULu2 = exp(log(var(cu))+z*se12)
  if(LLu2 < var & ULu2 > var) {
    smu2[i] = 1 } else {
      smu2[i] = 0}
  LLu3 = exp(log(var(cu))-z*se13)
  ULu3 = exp(log(var(cu))+z*se13)
  if(LLu3 < var & ULu3 > var) {
    smu3[i] = 1 } else {
      smu3[i] = 0}
} 
outu0 = mean(smu0)
outu1 = mean(smu1)
outu2 = mean(smu2)
outu3 = mean(smu3)

# CHI-QUADRO(1)
g = 1
var = 2 * g
smc0 = matrix(NA,M)
smc1 = matrix(NA,M)
smc2 = matrix(NA,M)
smc3 = matrix(NA,M)
for(i in 1:M) {
  cc = rchisq(n, g)
  perc = 1/(2*(n-4)^(1/2))
  gamma41 = n*sum((cc-mean(cc))^4)/(sum((cc-mean(cc))^2))^2
  gamma42 = n*sum((cc-mean(cc, perc, na.rm = F))^4)/(sum((cc-mean(cc))^2))^2
  gamma43 = n*sum((cc-median(cc))^4)/(sum((cc-mean(cc))^2))^2
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  se11 = ((gamma41-(n-3)/n)/(n-1))^(1/2)
  se12 = ((gamma42-(n-3)/n)/(n-1))^(1/2)
  se13 = ((gamma43-(n-3)/n)/(n-1))^(1/2)
  LLc0 = (n-1)*var(cc)/qchisq((1-alpha/2), n-1)
  ULc0 = (n-1)*var(cc)/qchisq((alpha/2), n-1)
  if(LLc0 < var & ULc0 > var) {
    smc0[i] = 1 } else {
      smc0[i] = 0}
  LLc1 = exp(log(var(cc))-z*se11)
  ULc1 = exp(log(var(cc))+z*se11)
  if(LLc1 < var & ULc1 > var) {
    smc1[i] = 1 } else {
      smc1[i] = 0}
  LLc2 = exp(log(var(cc))-z*se12)
  ULc2 = exp(log(var(cc))+z*se12)
  if(LLc2 < var & ULc2 > var) {
    smc2[i] = 1 } else {
      smc2[i] = 0}
  LLc3 = exp(log(var(cc))-z*se13)
  ULc3 = exp(log(var(cc))+z*se13)
  if(LLc3 < var & ULc3 > var) {
    smc3[i] = 1 } else {
      smc3[i] = 0}
} 
outc0 = mean(smc0)
outc1 = mean(smc1)
outc2 = mean(smc2)
outc3 = mean(smc3)

# T(5)
n = 10
M = 10000
g =  5
var = g/(g-2) 
smt0 = matrix(NA,M)
smt1 = matrix(NA,M)
smt2 = matrix(NA,M)
smt3 = matrix(NA,M)
for(i in 1:M) {
  ct = rt(n, g)
  perc = 1/(2*(n-4)^(1/2))
  gamma41 = n*sum((ct-mean(ct))^4)/(sum((ct-mean(ct))^2))^2
  gamma42 = n*sum((ct-mean(ct, perc, na.rm = F))^4)/(sum((ct-mean(ct))^2))^2
  gamma43 = n*sum((ct-median(ct))^4)/(sum((ct-mean(ct))^2))^2
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  se11 = ((gamma41-(n-3)/n)/(n-1))^(1/2)
  se12 = ((gamma42-(n-3)/n)/(n-1))^(1/2)
  se13 = ((gamma43-(n-3)/n)/(n-1))^(1/2)
  LLt0 = (n-1)*var(ct)/qchisq((1-alpha/2), n-1)
  ULt0 = (n-1)*var(ct)/qchisq((alpha/2), n-1)
  if(LLt0 < var & ULt0 > var) {
    smt0[i] = 1 } else {
      smt0[i] = 0}
  LLt1 = exp(log(var(ct))-z*se11)
  ULt1 = exp(log(var(ct))+z*se11)
  if(LLt1 < var & ULt1 > var) {
    smt1[i] = 1 } else {
      smt1[i] = 0}
  LLt2 = exp(log(var(ct))-z*se12)
  ULt2 = exp(log(var(ct))+z*se12)
  if(LLt2 < var & ULt2 > var) {
    smt2[i] = 1 } else {
      smt2[i] = 0}
  LLt3 = exp(log(var(ct))-z*se13)
  ULt3 = exp(log(var(ct))+z*se13)
  if(LLt3 < var & ULt3 > var) {
    smt3[i] = 1 } else {
      smt3[i] = 0}
} 
outt0 = mean(smt0)
outt1 = mean(smt1)
outt2 = mean(smt2)
outt3 = mean(smt3)

# ESPONENZIALE(1)
n = 10
M = 10000
var = 1
sme0 = matrix(NA,M)
sme1 = matrix(NA,M)
sme2 = matrix(NA,M)
sme3 = matrix(NA,M)
for(i in 1:M) {
  ce = rexp(n)
  perc = 1/(2*(n-4)^(1/2))
  gamma41 = n*sum((ce-mean(ce))^4)/(sum((ce-mean(ce))^2))^2
  gamma42 = n*sum((ce-mean(ce, perc, na.rm = F))^4)/(sum((ce-mean(ce))^2))^2
  gamma43 = n*sum((ce-median(ce))^4)/(sum((ce-mean(ce))^2))^2
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  se11 = ((gamma41-(n-3)/n)/(n-1))^(1/2)
  se12 = ((gamma42-(n-3)/n)/(n-1))^(1/2)
  se13 = ((gamma43-(n-3)/n)/(n-1))^(1/2)
  LLe0 = (n-1)*var(ce)/qchisq((1-alpha/2), n-1)
  ULe0 = (n-1)*var(ce)/qchisq((alpha/2), n-1)
  if(LLe0 < var & ULe0 > var) {
    sme0[i] = 1 } else {
      sme0[i] = 0}
  LLe1 = exp(log(var(ce))-z*se11)
  ULe1 = exp(log(var(ce))+z*se11)
  if(LLe1 < var & ULe1 > var) {
    sme1[i] = 1 } else {
      sme1[i] = 0}
  LLe2 = exp(log(var(ce))-z*se12)
  ULe2 = exp(log(var(ce))+z*se12)
  if(LLe2 < var & ULe2 > var) {
    sme2[i] = 1 } else {
      sme2[i] = 0}
  LLe3 = exp(log(var(ce))-z*se13)
  ULe3 = exp(log(var(ce))+z*se13)
  if(LLe3 < var & ULe3 > var) {
    sme3[i] = 1 } else {
      sme3[i] = 0}
} 
oute0 = mean(sme0)
oute1 = mean(sme1)
oute2 = mean(sme2)
oute3 = mean(sme3)

# BETA(3,3)
a = 3
b = 3
var = (a*b)/((a+b)^2*(a+b+1))
sm0b = matrix(NA, M)
sm1b = matrix(NA, M)
sm2b = matrix(NA, M)
sm3b = matrix(NA, M)
for(i in 1:M){
  camp = rbeta(n, a, b)
  perc = 1/(2*(n-4)^(1/2))
  gamma41 = n*sum((camp-mean(camp))^4)/(sum((camp-mean(camp))^2))^2
  gamma42 = n*sum((camp-mean(camp, perc, na.rm = F))^4)/(sum((camp-mean(camp))^2))^2
  gamma43 = n*sum((camp-median(camp))^4)/(sum((camp-mean(camp))^2))^2
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  se11 = ((gamma41-(n-3)/n)/(n-1))^(1/2)
  se12 = ((gamma42-(n-3)/n)/(n-1))^(1/2)
  se13 = ((gamma43-(n-3)/n)/(n-1))^(1/2)
  
  LL0 = (n-1)*var(camp)/qchisq((1-alpha/2), n-1)
  UL0 = (n-1)*var(camp)/qchisq((alpha/2), n-1)
  if(LL0 < var & UL0 > var) {
    sm0b[i] = 1 } else {
      sm0b[i] = 0}
  LL1 = exp(log(var(camp))-z*se11)
  UL1 = exp(log(var(camp)) + z * se11)
  if(LL1 < var & UL1 > var) {
    sm1b[i] = 1 } else{
      sm1b[i] = 0}
  LL2 = exp(log(var(camp)) - z * se12)
  UL2 = exp(log(var(camp)) + z * se12)
  if(LL2 < var & UL2 > var) {
    sm2b[i] = 1} else {
      sm2b[i] = 0}
  LL3 = exp(log(var(camp)) - z * se13)
  UL3 = exp(log(var(camp)) + z * se13)
  if(LL3 < var & UL3 > var) {
    sm3b[i] = 1} else {
      sm3b[i] = 0}
}
out0 = mean(sm0b)
out1 = mean(sm1b)
out2 = mean(sm2b)
out3 = mean(sm3b)

################################################################################
# TABLE 4: STIMA PROBABILITA DI COPERTURA BASATA SUGLI STIMAOTRI DELLA CURTOSI CON SE(2)

# N(0,1)
sm1 = matrix(NA, M)
sm2 = matrix(NA, M)
sm3 = matrix(NA, M)
for(i in 1:M){
  camp = rnorm(n, 0, 1)
  var = 1
  perc = 1/(2*(n-4)^(1/2))
  gamma41 = n*sum((camp-mean(camp))^4)/(sum((camp-mean(camp))^2))^2
  gamma42 = n*sum((camp-mean(camp, perc/2, na.rm = F))^4)/(sum((camp-mean(camp))^2))^2
  gamma43 = n*sum((camp-median(camp))^4)/(sum((camp-mean(camp))^2))^2
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c = n/(n-z)
  se21 = c*((gamma41-(n-3)/n)/(n-1))^(1/2)
  se22 = c*((gamma42-(n-3)/n)/(n-1))^(1/2)
  se23 = c*((gamma43-(n-3)/n)/(n-1))^(1/2)
  LL1 = exp(log(var(camp))-z*se21)
  UL1 = exp(log(var(camp)) + z * se21)
  if(LL1 < var & UL1 > var) {
    sm1[i] = 1 } else{
      sm1[i] = 0}
  LL2 = exp(log(var(camp)) - z * se22)
  UL2 = exp(log(var(camp)) + z * se22)
  if(LL2 < var & UL2 > var) {
    sm2[i] = 1} else {
      sm2[i] = 0}
  LL3 = exp(log(var(camp)) - z * se23)
  UL3 = exp(log(var(camp)) + z * se23)
  if(LL3 < var & UL3 > var) {
    sm3[i] = 1} else {
      sm3[i] = 0}
}
out1 = mean(sm1)
out2 = mean(sm2)
out3 = mean(sm3)

# U(0,1)
smu1 = matrix(NA,M)
smu2 = matrix(NA,M)
smu3 = matrix(NA,M)
for(i in 1:M) {
  cu = runif(n, 0, 1)
  var = 1/12
  perc = 1/(2*(n-4)^(1/2))
  gamma41 = n*sum((cu-mean(cu))^4)/(sum((cu-mean(cu))^2))^2
  gamma42 = n*sum((cu-mean(cu, perc, na.rm = F))^4)/(sum((cu-mean(cu))^2))^2
  gamma43 = n*sum((cu-median(cu))^4)/(sum((cu-mean(cu))^2))^2
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c = n/(n-z)
  se21 = c*((gamma41-(n-3)/n)/(n-1))^(1/2)
  se22 = c*((gamma42-(n-3)/n)/(n-1))^(1/2)
  se23 = c*((gamma43-(n-3)/n)/(n-1))^(1/2)
  LLu1 = exp(log(var(cu))-z*se21)
  ULu1 = exp(log(var(cu))+z*se21)
  if(LLu1 < var & ULu1 > var) {
    smu1[i] = 1 } else {
      smu1[i] = 0}
  LLu2 = exp(log(var(cu))-z*se22)
  ULu2 = exp(log(var(cu))+z*se22)
  if(LLu2 < var & ULu2 > var) {
    smu2[i] = 1 } else {
      smu2[i] = 0}
  LLu3 = exp(log(var(cu))-z*se23)
  ULu3 = exp(log(var(cu))+z*se23)
  if(LLu3 < var & ULu3 > var) {
    smu3[i] = 1 } else {
      smu3[i] = 0}
} 
outu1 = mean(smu1)
outu2 = mean(smu2)
outu3 = mean(smu3)

# CHI-QUADRO(1)
g = 1
var = 2 * g
smc1 = matrix(NA,M)
smc2 = matrix(NA,M)
smc3 = matrix(NA,M)
for(i in 1:M) {
  cc = rchisq(n, g)
  perc = 1/(2*(n-4)^(1/2))
  gamma41 = n*sum((cc-mean(cc))^4)/(sum((cc-mean(cc))^2))^2
  gamma42 = n*sum((cc-mean(cc, perc, na.rm = F))^4)/(sum((cc-mean(cc))^2))^2
  gamma43 = n*sum((cc-median(cc))^4)/(sum((cc-mean(cc))^2))^2
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c = n/(n-z)
  se21 = c*((gamma41-(n-3)/n)/(n-1))^(1/2)
  se22 = c*((gamma42-(n-3)/n)/(n-1))^(1/2)
  se23 = c*((gamma43-(n-3)/n)/(n-1))^(1/2)
  LLc1 = exp(log(var(cc))-z*se21)
  ULc1 = exp(log(var(cc))+z*se21)
  if(LLc1 < var & ULc1 > var) {
    smc1[i] = 1 } else {
      smc1[i] = 0}
  LLc2 = exp(log(var(cc))-z*se22)
  ULc2 = exp(log(var(cc))+z*se22)
  if(LLc2 < var & ULc2 > var) {
    smc2[i] = 1 } else {
      smc2[i] = 0}
  LLc3 = exp(log(var(cc))-z*se23)
  ULc3 = exp(log(var(cc))+z*se23)
  if(LLc3 < var & ULc3 > var) {
    smc3[i] = 1 } else {
      smc3[i] = 0}
} 
outc1 = mean(smc1)
outc2 = mean(smc2)
outc3 = mean(smc3)

# T-STUDENT(5)
n = 10
M = 10000
g = 5
var = g/(g-2)
smt1 = matrix(NA,M)
smt2 = matrix(NA,M)
smt3 = matrix(NA,M)
for(i in 1:M) {
  ct = rt(n, g)
  perc = 1/(2*(n-4)^(1/2))
  gamma41 = n*sum((ct-mean(ct))^4)/(sum((ct-mean(ct))^2))^2
  gamma42 = n*sum((ct-mean(ct, perc, na.rm = F))^4)/(sum((ct-mean(ct))^2))^2
  gamma43 = n*sum((ct-median(ct))^4)/(sum((ct-mean(ct))^2))^2
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c = n/(n-z)
  se21 = c*((gamma41-(n-3)/n)/(n-1))^(1/2)
  se22 = c*((gamma42-(n-3)/n)/(n-1))^(1/2)
  se23 = c*((gamma43-(n-3)/n)/(n-1))^(1/2)
  LLt1 = exp(log(var(ct))-z*se21)
  ULt1 = exp(log(var(ct))+z*se21)
  if(LLt1 < var & ULt1 > var) {
    smt1[i] = 1 } else {
      smt1[i] = 0}
  LLt2 = exp(log(var(ct))-z*se22)
  ULt2 = exp(log(var(ct))+z*se22)
  if(LLt2 < var & ULt2 > var) {
    smt2[i] = 1 } else {
      smt2[i] = 0}
  LLt3 = exp(log(var(ct))-z*se23)
  ULt3 = exp(log(var(ct))+z*se23)
  if(LLt3 < var & ULt3 > var) {
    smt3[i] = 1 } else {
      smt3[i] = 0}
} 
outt1 = mean(smt1)
outt2 = mean(smt2)
outt3 = mean(smt3)

# ESPONENZIALE(1)
var = 1
sme1 = matrix(NA,M)
sme2 = matrix(NA,M)
sme3 = matrix(NA,M)
for(i in 1:M) {
  ce = rexp(n)
  perc = 1/(2*(n-4)^(1/2))
  gamma41 = n*sum((ce-mean(ce))^4)/(sum((ce-mean(ce))^2))^2
  gamma42 = n*sum((ce-mean(ce, perc, na.rm = F))^4)/(sum((ce-mean(ce))^2))^2
  gamma43 = n*sum((ce-median(ce))^4)/(sum((ce-mean(ce))^2))^2
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c = n/(n-z)
  se21 = c*((gamma41-(n-3)/n)/(n-1))^(1/2)
  se22 = c*((gamma42-(n-3)/n)/(n-1))^(1/2)
  se23 = c*((gamma43-(n-3)/n)/(n-1))^(1/2)
  LLe1 = exp(log(var(ce))-z*se21)
  ULe1 = exp(log(var(ce))+z*se21)
  if(LLe1 < var & ULe1 > var) {
    sme1[i] = 1 } else {
      sme1[i] = 0}
  LLe2 = exp(log(var(ce))-z*se22)
  ULe2 = exp(log(var(ce))+z*se22)
  if(LLe2 < var & ULe2 > var) {
    sme2[i] = 1 } else {
      sme2[i] = 0}
  LLe3 = exp(log(var(ce))-z*se23)
  ULe3 = exp(log(var(ce))+z*se23)
  if(LLe3 < var & ULe3 > var) {
    sme3[i] = 1 } else {
      sme3[i] = 0}
} 
oute1 = mean(sme1)
oute2 = mean(sme2)
oute3 = mean(sme3)

# BETA(3,3)
a = 3
b = 3
var = (a*b)/((a+b)^2*(a+b+1))
sm1b = matrix(NA, M)
sm2b = matrix(NA, M)
sm3b = matrix(NA, M)
for(i in 1:M){
  camp = rbeta(n, a, b)
  perc = 1/(2*(n-4)^(1/2))
  gamma41 = n*sum((camp-mean(camp))^4)/(sum((camp-mean(camp))^2))^2
  gamma42 = n*sum((camp-mean(camp, perc/2, na.rm = F))^4)/(sum((camp-mean(camp))^2))^2
  gamma43 = n*sum((camp-median(camp))^4)/(sum((camp-mean(camp))^2))^2
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c = n/(n-z)
  se21 = c*((gamma41-(n-3)/n)/(n-1))^(1/2)
  se22 = c*((gamma42-(n-3)/n)/(n-1))^(1/2)
  se23 = c*((gamma43-(n-3)/n)/(n-1))^(1/2)
  LL1 = exp(log(var(camp))-z*se21)
  UL1 = exp(log(var(camp)) + z * se21)
  if(LL1 < var & UL1 > var) {
    sm1b[i] = 1 } else{
      sm1b[i] = 0}
  LL2 = exp(log(var(camp)) - z * se22)
  UL2 = exp(log(var(camp)) + z * se22)
  if(LL2 < var & UL2 > var) {
    sm2b[i] = 1} else {
      sm2b[i] = 0}
  LL3 = exp(log(var(camp)) - z * se23)
  UL3 = exp(log(var(camp)) + z * se23)
  if(LL3 < var & UL3 > var) {
    sm3b[i] = 1} else {
      sm3b[i] = 0}
}
out1 = mean(sm1b)
out2 = mean(sm2b)
out3 = mean(sm3b)


################################################################################
# RAPPORTO DI VARIANZE
################################################################################
# TABELLA 5: STIMA PROBABILITA DI COPERTURA PER IL RAPPORTO DI VARIANZE

# N(0,1)
var = 1
n = 10
M = 10000
smf1 = matrix(NA, M)
smn = matrix(NA, M)
sm1 = matrix(NA, M)
sm2 = matrix(NA, M)
sm3 = matrix(NA, M)
for(i in 1:M){
  c1 = rnorm(n, 0, 1)
  c2 = rnorm(n, 0, 1)
  perc = 1/(2*(n-4)^(1/2))
  gamma41_1 = n*sum((c1-mean(c1))^4)/(sum((c1-mean(c1))^2))^2
  gamma41_2 = n*sum((c2-mean(c2))^4)/(sum((c2-mean(c2))^2))^2
  gamma42_1 = n*sum((c1-mean(c1, perc/2, na.rm = F))^4)/(sum((c1-mean(c1))^2))^2
  gamma42_2 = n*sum((c2-mean(c2, perc/2, na.rm = F))^4)/(sum((c2-mean(c2))^2))^2
  gamma43_1 = n*sum((c1-median(c1))^4)/(sum((c1-mean(c1))^2))^2
  gamma43_2 = n*sum((c2-median(c2))^4)/(sum((c2-mean(c2))^2))^2
  
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c = n/(n-z)
  # se di table 3
  se11_1 = c * ((gamma41_1 - (n-3)/n)/(n-1))^(1/2)
  se11_2 = c * ((gamma41_2 - (n-3)/n)/(n-1))^(1/2)
  se12_1 = c * ((gamma42_1 - (n-3)/n)/(n-1))^(1/2)
  se12_2 = c * ((gamma42_2 - (n-3)/n)/(n-1))^(1/2)
  se13_1 = c * ((gamma43_1 - (n-3)/n)/(n-1))^(1/2)
  se13_2 = c * ((gamma43_2 - (n-3)/n)/(n-1))^(1/2)
  # normal
  LLn = (var(c2)/var(c1)) * qf((alpha/2), n-1, n-1)
  ULn = (var(c2)/var(c1)) * qf((1-alpha/2), n-1, n-1)
  if(LLn < var & ULn > var) {
    smn[i] = 1} else{
      smn[i] = 0}
  # F1
  f11 = sum((c1 - mean(c1))^4)
  f12 = sum((c2 - mean(c2))^4)
  mu4f1 = (f11 + f12)/(n+n)
  s2f1 = ((n-1) * var(c1) + (n-1) * var(c2))/(n+n)
  gamma4f1 = mu4f1/(s2f1^2)
  r1 = (2*n) / (gamma4f1 - ((n-3)/(n-1)))
  r2 = (2*n) / (gamma4f1 - ((n-3)/(n-1)))
  # intervallo F1
  LLf1 = (var(c2)/var(c1)) * qf((alpha/2), r1, r2)
  ULf1 = (var(c2)/var(c1)) * qf((1-alpha/2), r1, r2)
  if(LLf1 < var & ULf1 > var) {
    smf1[i] = 1} else {
      smf1[i] = 0}
  
  LL1 = exp(log(c*var(c2))-log(c*var(c1))-z*(se11_1+se11_2))
  UL1 = exp(log(c*var(c2))-log(c*var(c1))+z*(se11_1+se11_2))
  if(LL1 < var & UL1 > var) {
  sm1[i] = 1 } else {
      sm1[i] = 0}
  LL2 = exp(log(c*var(c2))-log(c*var(c1))-z*(se12_1+se12_2))
  UL2 = exp(log(c*var(c2))-log(c*var(c1))+z*(se12_1+se12_2))
  if(LL2 < var & UL2 > var) {
    sm2[i] = 1 } else {
      sm2[i] = 0}
  LL3 = exp(log(c*var(c2))-log(c*var(c1))-z*(se13_1+se13_2))
  UL3 = exp(log(c*var(c2))-log(c*var(c1))+z*(se13_1+se13_2))
  if(LL3 < var & UL3 > var) {
    sm3[i] = 1 } else {
      sm3[i] = 0}
}
outF1 = mean(smf1)
outN = mean(smn)
out1 = mean(sm1)
out2 = mean(sm2)
out3 = mean(sm3)

# U(0,1)
var = 1
smf1u = matrix(NA, M)
smnu = matrix(NA, M)
sm1u = matrix(NA, M)
sm2u = matrix(NA, M)
sm3u = matrix(NA, M)
for(i in 1:M){
  c1 = runif(n, 0, 1)
  c2 = runif(n, 0, 1)
  perc = 1/(2*(n-4)^(1/2))
  
  gamma41_1 = n*sum((c1-mean(c1))^4)/(sum((c1-mean(c1))^2))^2
  gamma41_2 = n*sum((c2-mean(c2))^4)/(sum((c2-mean(c2))^2))^2
  gamma42_1 = n*sum((c1-mean(c1, perc/2, na.rm = F))^4)/(sum((c1-mean(c1))^2))^2
  gamma42_2 = n*sum((c2-mean(c2, perc/2, na.rm = F))^4)/(sum((c2-mean(c2))^2))^2
  gamma43_1 = n*sum((c1-median(c1))^4)/(sum((c1-mean(c1))^2))^2
  gamma43_2 = n*sum((c2-median(c2))^4)/(sum((c2-mean(c2))^2))^2
  
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c = n/(n-z)
  # se di table 3
  se11_1 = c * ((gamma41_1 - (n-3)/n)/(n-1))^(1/2)
  se11_2 = c * ((gamma41_2 - (n-3)/n)/(n-1))^(1/2)
  se12_1 = c * ((gamma42_1 - (n-3)/n)/(n-1))^(1/2)
  se12_2 = c * ((gamma42_2 - (n-3)/n)/(n-1))^(1/2)
  se13_1 = c * ((gamma43_1 - (n-3)/n)/(n-1))^(1/2)
  se13_2 = c * ((gamma43_2 - (n-3)/n)/(n-1))^(1/2)
  # normal
  LLnu = (var(c2)/var(c1)) * qf((alpha/2), n-1, n-1)
  ULnu = (var(c2)/var(c1)) * qf((1-alpha/2), n-1, n-1)
  if(LLnu < var & ULnu > var) {
    smnu[i] = 1} else{
      smnu[i] = 0}
  # F1
  f11 = sum((c1 - mean(c1))^4)
  f12 = sum((c2 - mean(c2))^4)
  mu4f1 = (f11 + f12)/(2*n)
  s2f1 = ((n-1) * var(c1) + (n-1) * var(c2))/(n+n)
  gamma4f1 = mu4f1/(s2f1^2)
  r1 = (2*n) / (gamma4f1 - ((n-3)/(n-1)))
  r2 = (2*n) / (gamma4f1 - ((n-3)/(n-1)))
  # intervallo F1
  LLf1u = (var(c2)/var(c1)) * qf((alpha/2), r1, r2)
  ULf1u = (var(c2)/var(c1)) * qf((1-alpha/2), r1, r2)
  if(LLf1u < var & ULf1u > var) {
    smf1u[i] = 1} else {
      smf1u[i] = 0}
  #gamma4
  LL1u = exp(log(c*var(c2))-log(c*var(c1))-z*(se11_1+se11_2))
  UL1u = exp(log(c*var(c2))-log(c*var(c1))+z*(se11_1+se11_2))
  if(LL1u < var & UL1u > var) {
    sm1u[i] = 1 } else {
      sm1u[i] = 0}
  LL2u = exp(log(c*var(c2))-log(c*var(c1))-z*(se12_1+se12_2))
  UL2u = exp(log(c*var(c2))-log(c*var(c1))+z*(se12_1+se12_2))
  if(LL2u < var & UL2u > var) {
    sm2u[i] = 1 } else {
      sm2u[i] = 0}
  LL3u = exp(log(c*var(c2))-log(c*var(c1))-z*(se13_1+se13_2))
  UL3u = exp(log(c*var(c2))-log(c*var(c1))+z*(se13_1+se13_2))
  if(LL3u < var & UL3u > var) {
    sm3u[i] = 1 } else {
      sm3u[i] = 0}
}
outF1u = mean(smf1u)
outNu = mean(smnu)
out1u = mean(sm1u)
out2u = mean(sm2u)
out3u = mean(sm3u)

# CHI-QUADRO(1)
g = 1
var = 2 * g
smfc = matrix(NA, M)
smnc = matrix(NA, M)
sm1c = matrix(NA, M)
sm2c = matrix(NA, M)
sm3c = matrix(NA, M)
for(i in 1:M){
  c1 = rchisq(n, g)
  c2 = rchisq(n, g)
  tp = 1/(2*(n-4)^(1/2))
  gamma41_1 = n*sum((c1-mean(c1))^4)/(sum((c1-mean(c1))^2))^2
  gamma41_2 = n*sum((c2-mean(c2))^4)/(sum((c2-mean(c2))^2))^2
  gamma42_1 = n*sum((c1-mean(c1, tp/2, na.rm = F))^4)/(sum((c1-mean(c1))^2))^2
  gamma42_2 = n*sum((c2-mean(c2, tp/2, na.rm = F))^4)/(sum((c2-mean(c2))^2))^2
  gamma43_1 = n*sum((c1-median(c1))^4)/(sum((c1-mean(c1))^2))^2
  gamma43_2 = n*sum((c2-median(c2))^4)/(sum((c2-mean(c2))^2))^2
  
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c = n/(n-z)
  se11_1 = c * ((gamma41_1 - (n-3)/n)/(n-1))^(1/2)
  se11_2 = c * ((gamma41_2 - (n-3)/n)/(n-1))^(1/2)
  se12_1 = c * ((gamma42_1 - (n-3)/n)/(n-1))^(1/2)
  se12_2 = c * ((gamma42_2 - (n-3)/n)/(n-1))^(1/2)
  se13_1 = c * ((gamma43_1 - (n-3)/n)/(n-1))^(1/2)
  se13_2 = c * ((gamma43_2 - (n-3)/n)/(n-1))^(1/2)
  LLnc = (var(c2)/var(c1)) * qf((alpha/2), n-1, n-1)
  ULnc = (var(c2)/var(c1)) * qf((1-alpha/2), n-1, n-1)
  if(LLnc < var & ULnc > var) {
    smnc[i] = 1} else{
      smnc[i] = 0}
  f11 = sum((c1 - mean(c1))^4)
  f12 = sum((c2 - mean(c2))^4)
  mu4f1 = (f11 + f12)/(2*n)
  s2f1 = ((n-1) * var(c1) + (n-1) * var(c2))/(n+n)
  gamma4f1 = mu4f1/(s2f1^2)
  r1 = (2*n) / (gamma4f1 - ((n-3)/(n-1)))
  r2 = (2*n) / (gamma4f1 - ((n-3)/(n-1)))
  LLf1c = (var(c2)/var(c1)) * qf((alpha/2), r1, r2)
  ULf1c = (var(c2)/var(c1)) * qf((1-alpha/2), r1, r2)
  if(LLf1c < var & ULf1c > var) {
    smfc[i] = 1} else {
      smfc[i] = 0}
  
  LL1c = exp(log(c*var(c2))-log(c*var(c1))-z*(se11_1+se11_2))
  UL1c = exp(log(c*var(c2))-log(c*var(c1))+z*(se11_1+se11_2))
  if(LL1c < var & UL1c > var) {
    sm1c[i] = 1 } else {
      sm1c[i] = 0}
  LL2c = exp(log(c*var(c2))-log(c*var(c1))-z*(se12_1+se12_2))
  UL2c = exp(log(c*var(c2))-log(c*var(c1))+z*(se12_1+se12_2))
  if(LL2c < var & UL2c > var) {
    sm2c[i] = 1 } else {
      sm2c[i] = 0}
  LL3c = exp(log(c*var(c2))-log(c*var(c1))-z*(se13_1+se13_2))
  UL3c = exp(log(c*var(c2))-log(c*var(c1))+z*(se13_1+se13_2))
  if(LL3c < var & UL3c > var) {
    sm3c[i] = 1 } else {
      sm3c[i] = 0}
}
outF1c = mean(smfc)
outNc = mean(smnc)
out1c = mean(sm1c)
out2c = mean(sm2c)
out3c = mean(sm3c)

# T-STUDENT(5)
g = 5
var = 1
smft = matrix(NA, M)
smnt = matrix(NA, M)
sm1t = matrix(NA, M)
sm2t = matrix(NA, M)
sm3t = matrix(NA, M)
for(i in 1:M){
  c1 = rt(n, g)
  c2 = rt(n, g)
  perc = 1/(2*(n-4)^(1/2))
  gamma41_1 = n*sum((c1-mean(c1))^4)/(sum((c1-mean(c1))^2))^2
  gamma41_2 = n*sum((c2-mean(c2))^4)/(sum((c2-mean(c2))^2))^2
  gamma42_1 = n*sum((c1-mean(c1, perc/2, na.rm = F))^4)/(sum((c1-mean(c1))^2))^2
  gamma42_2 = n*sum((c2-mean(c2, perc/2, na.rm = F))^4)/(sum((c2-mean(c2))^2))^2
  gamma43_1 = n*sum((c1-median(c1))^4)/(sum((c1-mean(c1))^2))^2
  gamma43_2 = n*sum((c2-median(c2))^4)/(sum((c2-mean(c2))^2))^2
  
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c = n/(n-z)
  se11_1 = c * ((gamma41_1 - (n-3)/n)/(n-1))^(1/2)
  se11_2 = c * ((gamma41_2 - (n-3)/n)/(n-1))^(1/2)
  se12_1 = c * ((gamma42_1 - (n-3)/n)/(n-1))^(1/2)
  se12_2 = c * ((gamma42_2 - (n-3)/n)/(n-1))^(1/2)
  se13_1 = c * ((gamma43_1 - (n-3)/n)/(n-1))^(1/2)
  se13_2 = c * ((gamma43_2 - (n-3)/n)/(n-1))^(1/2)
  
  LLnt = (var(c2)/var(c1)) * qf((alpha/2), n-1, n-1)
  ULnt = (var(c2)/var(c1)) * qf((1-alpha/2), n-1, n-1)
  if(LLnt < var & ULnt > var) {
    smnt[i] = 1} else{
      smnt[i] = 0}
  f11 = sum((c1 - mean(c1))^4)
  f12 = sum((c2 - mean(c2))^4)
  mu4f1 = (f11 + f12)/(2*n)
  s2f1 = ((n-1) * var(c1) + (n-1) * var(c2))/(n+n)
  gamma4f1 = mu4f1/(s2f1^2)
  r1 = (2*n) / (gamma4f1 - ((n-3)/(n-1)))
  r2 = (2*n) / (gamma4f1 - ((n-3)/(n-1)))
  LLf1t = (var(c2)/var(c1)) * qf((alpha/2), r1, r2)
  ULf1t = (var(c2)/var(c1)) * qf((1-alpha/2), r1, r2)
  if(LLf1t < var & ULf1t > var) {
    smft[i] = 1} else {
      smft[i] = 0}
  
  LL1t = exp(log(c*var(c2))-log(c*var(c1))-z*(se11_1+se11_2))
  UL1t = exp(log(c*var(c2))-log(c*var(c1))+z*(se11_1+se11_2))
  if(LL1t < var & UL1t > var) {
    sm1t[i] = 1 } else {
      sm1t[i] = 0}
  LL2t = exp(log(c*var(c2))-log(c*var(c1))-z*(se12_1+se12_2))
  UL2t = exp(log(c*var(c2))-log(c*var(c1))+z*(se12_1+se12_2))
  if(LL2t < var & UL2t > var) {
    sm2t[i] = 1 } else {
      sm2t[i] = 0}
  LL3t = exp(log(c*var(c2))-log(c*var(c1))-z*(se13_1+se13_2))
  UL3t = exp(log(c*var(c2))-log(c*var(c1))+z*(se13_1+se13_2))
  if(LL3t < var & UL3t > var) {
    sm3t[i] = 1 } else {
      sm3t[i] = 0}
}
outF1t = mean(smft)
outNt = mean(smnt)
out1t = mean(sm1t)
out2t = mean(sm2t)
out3t = mean(sm3t)

# ESPONENZIALE(1)
var = 1
smfe = matrix(NA, M)
smne = matrix(NA, M)
sm1e = matrix(NA, M)
sm2e = matrix(NA, M)
sm3e = matrix(NA, M)
for(i in 1:M){
  c1 = rexp(n)
  c2 = rexp(n)
  perc = 1/(2*(n-4)^(1/2))
  gamma41_1 = n*sum((c1-mean(c1))^4)/(sum((c1-mean(c1))^2))^2
  gamma41_2 = n*sum((c2-mean(c2))^4)/(sum((c2-mean(c2))^2))^2
  gamma42_1 = n*sum((c1-mean(c1, perc/2, na.rm = F))^4)/(sum((c1-mean(c1))^2))^2
  gamma42_2 = n*sum((c2-mean(c2, perc/2, na.rm = F))^4)/(sum((c2-mean(c2))^2))^2
  gamma43_1 = n*sum((c1-median(c1))^4)/(sum((c1-mean(c1))^2))^2
  gamma43_2 = n*sum((c2-median(c2))^4)/(sum((c2-mean(c2))^2))^2
  
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c = n/(n-z)
  se11_1 = c * ((gamma41_1 - (n-3)/n)/(n-1))^(1/2)
  se11_2 = c * ((gamma41_2 - (n-3)/n)/(n-1))^(1/2)
  se12_1 = c * ((gamma42_1 - (n-3)/n)/(n-1))^(1/2)
  se12_2 = c * ((gamma42_2 - (n-3)/n)/(n-1))^(1/2)
  se13_1 = c * ((gamma43_1 - (n-3)/n)/(n-1))^(1/2)
  se13_2 = c * ((gamma43_2 - (n-3)/n)/(n-1))^(1/2)
  
  LLne = (var(c2)/var(c1)) * qf((alpha/2), n-1, n-1)
  ULne = (var(c2)/var(c1)) * qf(( 1-alpha/2), n-1, n-1)
  if(LLne < var & ULne > var) {
    smne[i] = 1} else{
      smne[i] = 0}
  f11 = sum((c1 - mean(c1))^4)
  f12 = sum((c2 - mean(c2))^4)
  mu4f1 = (f11 + f12)/(2*n)
  s2f1 = ((n-1) * var(c1) + (n-1) * var(c2))/(n+n)
  gamma4f1 = mu4f1/(s2f1^2)
  r1 = (2*n) / (gamma4f1 - ((n-3)/(n-1)))
  r2 = (2*n) / (gamma4f1 - ((n-3)/(n-1)))
  LLf1e = (var(c2)/var(c1)) * qf((alpha/2), r1, r2)
  ULf1e = (var(c2)/var(c1)) * qf((1-alpha/2), r1, r2)
  if(LLf1e < var & ULf1e > var) {
    smfe[i] = 1} else {
      smfe[i] = 0}
  
  LL1e = exp(log(c*var(c2))-log(c*var(c1))-z*(se11_1+se11_2))
  UL1e = exp(log(c*var(c2))-log(c*var(c1))+z*(se11_1+se11_2))
  if(LL1e < var & UL1e > var) {
    sm1e[i] = 1 } else {
      sm1e[i] = 0}
  LL2e = exp(log(c*var(c2))-log(c*var(c1))-z*(se12_1+se12_2))
  UL2e = exp(log(c*var(c2))-log(c*var(c1))+z*(se12_1+se12_2))
  if(LL2e < var & UL2e > var) {
    sm2e[i] = 1 } else {
      sm2e[i] = 0}
  LL3e = exp(log(c*var(c2))-log(c*var(c1))-z*(se13_1+se13_2))
  UL3e = exp(log(c*var(c2))-log(c*var(c1))+z*(se13_1+se13_2))
  if(LL3e < var & UL3e > var) {
    sm3e[i] = 1 } else {
      sm3e[i] = 0}
}
outF1e = mean(smfe)
outNe = mean(smne)
out1e = mean(sm1e)
out2e = mean(sm2e)
out3e = mean(sm3e)

#BETA(3,3)
a = 3
b = 3
var = 1
smf1b = matrix(NA, M)
smnb = matrix(NA, M)
sm1b = matrix(NA, M)
sm2b = matrix(NA, M)
sm3b = matrix(NA, M)
for(i in 1:M){
  c1 = rbeta(n, a, b)
  c2 = rbeta(n, a, b)
  perc = 1/(2*(n-4)^(1/2))
  gamma41_1 = n*sum((c1-mean(c1))^4)/(sum((c1-mean(c1))^2))^2
  gamma41_2 = n*sum((c2-mean(c2))^4)/(sum((c2-mean(c2))^2))^2
  gamma42_1 = n*sum((c1-mean(c1, perc/2, na.rm = F))^4)/(sum((c1-mean(c1))^2))^2
  gamma42_2 = n*sum((c2-mean(c2, perc/2, na.rm = F))^4)/(sum((c2-mean(c2))^2))^2
  gamma43_1 = n*sum((c1-median(c1))^4)/(sum((c1-mean(c1))^2))^2
  gamma43_2 = n*sum((c2-median(c2))^4)/(sum((c2-mean(c2))^2))^2
  
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c = n/(n-z)
  se11_1 = c * ((gamma41_1 - (n-3)/n)/(n-1))^(1/2)
  se11_2 = c * ((gamma41_2 - (n-3)/n)/(n-1))^(1/2)
  se12_1 = c * ((gamma42_1 - (n-3)/n)/(n-1))^(1/2)
  se12_2 = c * ((gamma42_2 - (n-3)/n)/(n-1))^(1/2)
  se13_1 = c * ((gamma43_1 - (n-3)/n)/(n-1))^(1/2)
  se13_2 = c * ((gamma43_2 - (n-3)/n)/(n-1))^(1/2)
  
  LLn = (var(c2)/var(c1)) * qf((alpha/2), n-1, n-1)
  ULn = (var(c2)/var(c1)) * qf((1-alpha/2), n-1, n-1)
  if(LLn < var & ULn > var) {
    smnb[i] = 1} else{
      smnb[i] = 0}
  f11 = sum((c1 - mean(c1))^4)
  f12 = sum((c2 - mean(c2))^4)
  mu4f1 = (f11 + f12)/(2*n)
  s2f1 = ((n-1) * var(c1) + (n-1) * var(c2))/(n+n)
  gamma4f1 = mu4f1/(s2f1^2)
  r1 = (2*n) / (gamma4f1 - ((n-3)/(n-1)))
  r2 = (2*n) / (gamma4f1 - ((n-3)/(n-1)))
  LLf1 = (var(c2)/var(c1)) * qf((alpha/2), r1, r2)
  ULf1 = (var(c2)/var(c1)) * qf((1-alpha/2), r1, r2)
  if(LLf1 < var & ULf1 > var) {
    smf1b[i] = 1} else {
      smf1b[i] = 0}
  LL1 = exp(log(c*var(c2))-log(c*var(c1))-z*(se11_1+se11_2))
  UL1 = exp(log(c*var(c2))-log(c*var(c1))+z*(se11_1+se11_2))
  if(LL1 < var & UL1 > var) {
    sm1b[i] = 1 } else {
      sm1b[i] = 0}
  LL2 = exp(log(c*var(c2))-log(c*var(c1))-z*(se12_1+se12_2))
  UL2 = exp(log(c*var(c2))-log(c*var(c1))+z*(se12_1+se12_2))
  if(LL2 < var & UL2 > var) {
    sm2b[i] = 1 } else {
      sm2b[i] = 0}
  LL3 = exp(log(c*var(c2))-log(c*var(c1))-z*(se13_1+se13_2))
  UL3 = exp(log(c*var(c2))-log(c*var(c1))+z*(se13_1+se13_2))
  if(LL3 < var & UL3 > var) {
    sm3b[i] = 1 } else {
      sm3b[i] = 0}
}
outF1 = mean(smf1b)
outN = mean(smnb)
out1 = mean(sm1b)
out2 = mean(sm2b)
out3 = mean(sm3b)


################################################################################
# MEDIA
################################################################################
# TABELLA 6: INTERVALLI DI CONFIDENZA PER LA MEDIA

#LN(0,1)
n = 10
M = 1000
mu = exp(1/2)
mug_x = function(x){
  out = exp(1/n * sum(log(x)))
  return(out)}
sim1 = matrix(NA, M)
sim2 = matrix(NA, M)
sim3 = matrix(NA, M)
sim4 = matrix(NA, M)
mux = matrix(NA, M)
mug = matrix(NA, M)
for(i in 1:M){
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c = n/(n-z)
  for(j in 1:M){
    camp = rlnorm(n, 0, 1)  
    mux[j] = mean(camp)
    mug[j] = mug_x(camp)
  }
  mod1 = lm(mux ~ 0 + mug)
  gamma = coef(mod1)
  # Normal
  t = qt(1-alpha/2, n-1)
  LL1 = mean(camp) - t * sd(camp)/(sqrt(n))
  UL1 = mean(camp) + t * sd(camp)/(sqrt(n))
  if(LL1 < mu & UL1 > mu) {
    sim1[i] = 1 } else {
      sim1[i] = 0}
  # Bonett
  sigmab = sqrt(exp(log(c * var(camp))))
  LL2 = mean(camp) - z * sigmab/(sqrt(n))
  UL2 = mean(camp) + z * sigmab/(sqrt(n))
  if(LL2 < mu & UL2 > mu) {
    sim2[i] = 1} else{
      sim2[i] = 0}
  # Johnson
  mu3 = mean((camp - mean(camp))^3) 
  LL3 = (mean(camp) + mu3/(6*var(camp)*n)) - t * sd(camp)/(sqrt(n))
  UL3 = (mean(camp) + mu3/(6*var(camp)*n)) + t * sd(camp)/(sqrt(n))
  if(LL3 < mu & UL3 > mu) {
    sim3[i] = 1} else{
      sim3[i] = 0}
  # GM
  LL4 = gamma * exp(mean(log(camp)) - t * sd(log(camp))/(sqrt(n))) 
  UL4 = gamma * exp(mean(log(camp)) + t * sd(log(camp))/(sqrt(n))) 
  if(LL4 < mu & UL4 > mu){
    sim4[i] = 1} else{
      sim4[i] = 0}
}
out1 = mean(sim1)
out2 = mean(sim2)
out3 = mean(sim3)
out4 = mean(sim4)

# U(0,1)
n = 10
M = 1000
mu = 1/2
mug_x = function(x){
  out = exp(1/n * sum(log(x)))
  return(out)}
sim1 = matrix(NA, M)
sim2 = matrix(NA, M)
sim3 = matrix(NA, M)
sim4 = matrix(NA, M)
mux = matrix(NA, M)
mug = matrix(NA, M)
for(i in 1:M){
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c = n/(n-z)
  for(j in 1:M){
    camp = runif(n, 0, 1)  
    mux[j] = mean(camp)
    mug[j] = mug_x(camp)
  }
  mod1 = lm(mux ~ 0 + mug)
  gamma = coef(mod1)
  # Normal
  t = qt(1-alpha/2, n-1)
  LL1 = mean(camp) - t * sd(camp)/(sqrt(n))
  UL1 = mean(camp) + t * sd(camp)/(sqrt(n))
  if(LL1 < mu & UL1 > mu) {
    sim1[i] = 1 } else {
      sim1[i] = 0}
  # Bonett
  sigmab = sqrt(exp(log(c * var(camp))))
  LL2 = mean(camp) - z * sigmab/(sqrt(n))
  UL2 = mean(camp) + z * sigmab/(sqrt(n))
  if(LL2 < mu & UL2 > mu) {
    sim2[i] = 1} else{
      sim2[i] = 0}
  # Johnson
  mu3 = mean((camp - mean(camp))^3) 
  LL3 = (mean(camp) + mu3/(6*var(camp)*n)) - t * sd(camp)/(sqrt(n))
  UL3 = (mean(camp) + mu3/(6*var(camp)*n)) + t * sd(camp)/(sqrt(n))
  if(LL3 < mu & UL3 > mu) {
    sim3[i] = 1} else{
      sim3[i] = 0}
  # GM
  LL4 = gamma * exp(mean(log(camp)) - t * sd(log(camp))/(sqrt(n))) 
  UL4 = gamma * exp(mean(log(camp)) + t * sd(log(camp))/(sqrt(n))) 
  if(LL4 < mu & UL4 > mu){
    sim4[i] = 1} else{
      sim4[i] = 0}
}
out1 = mean(sim1)
out2 = mean(sim2)
out3 = mean(sim3)
out4 = mean(sim4)

# EXP(1)
mu = 1
mug_x = function(x){
  out = exp(1/n * sum(log(x)))
  return(out)}
sim1 = matrix(NA, M)
sim2 = matrix(NA, M)
sim3 = matrix(NA, M)
sim4 = matrix(NA, M)
mux = matrix(NA, M)
mug = matrix(NA, M)
for(i in 1:M){
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c = n/(n-z)
  for(j in 1:M){
    camp = rexp(n,1)  
    mux[j] = mean(camp)
    mug[j] = mug_x(camp)
  }
  # Normal
  mod1 = lm(mux ~ 0 + mug)
  gamma = coef(mod1)
  t = qt(1-alpha/2, n-1)
  LL1 = mean(camp) - t * sd(camp)/(sqrt(n))
  UL1 = mean(camp) + t * sd(camp)/(sqrt(n))
  if(LL1 < mu & UL1 > mu) {
    sim1[i] = 1 } else {
      sim1[i] = 0}
  # Bonett
  sigmab = sqrt(exp(log(c * var(camp))))
  LL2 = mean(camp) - z * sigmab/(sqrt(n))
  UL2 = mean(camp) + z * sigmab/(sqrt(n))
  if(LL2 < mu & UL2 > mu) {
    sim2[i] = 1} else{
      sim2[i] = 0}
  # Johnson
  mu3 = mean((camp - mean(camp))^3) 
  LL3 = (mean(camp) + mu3/(6*var(camp)*n)) - t * sd(camp)/(sqrt(n))
  UL3 = (mean(camp) + mu3/(6*var(camp)*n)) + t * sd(camp)/(sqrt(n))
  if(LL3 < mu & UL3 > mu) {
    sim3[i] = 1} else{
      sim3[i] = 0}
  # GM
  LL4 = gamma * exp(mean(log(camp)) - t * sd(log(camp))/(sqrt(n))) 
  UL4 = gamma * exp(mean(log(camp)) + t * sd(log(camp))/(sqrt(n))) 
  if(LL4 < mu & UL4 > mu){
    sim4[i] = 1} else{
      sim4[i] = 0}
}
out1 = mean(sim1)
out2 = mean(sim2)
out3 = mean(sim3)
out4 = mean(sim4)

# CHISQ(1)
n = 10
M = 1000
mu = 1
mug_x = function(x){
  out = exp(1/n * sum(log(x)))
  return(out)}
sim1 = matrix(NA, M)
sim2 = matrix(NA, M)
sim3 = matrix(NA, M)
sim4 = matrix(NA, M)
mux = matrix(NA, M)
mug = matrix(NA, M)
for(i in 1:M){
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c = n/(n-z)
  for(j in 1:M){
    camp = rchisq(n, 1)  
    mux[j] = mean(camp)
    mug[j] = mug_x(camp)
  }
  # Normal
  mod1 = lm(mux ~ 0 + mug)
  gamma = coef(mod1)
  t = qt(1-alpha/2, n-1)
  LL1 = mean(camp) - t * sd(camp)/(sqrt(n))
  UL1 = mean(camp) + t * sd(camp)/(sqrt(n))
  if(LL1 < mu & UL1 > mu) {
    sim1[i] = 1 } else {
      sim1[i] = 0}
  # Bonett
  sigmab = sqrt(exp(log(c * var(camp))))
  LL2 = mean(camp) - z * sigmab/(sqrt(n))
  UL2 = mean(camp) + z * sigmab/(sqrt(n))
  if(LL2 < mu & UL2 > mu) {
    sim2[i] = 1} else{
      sim2[i] = 0}
  # Johnson
  mu3 = mean((camp - mean(camp))^3) 
  LL3 = (mean(camp) + mu3/(6*var(camp)*n)) - t * sd(camp)/(sqrt(n))
  UL3 = (mean(camp) + mu3/(6*var(camp)*n)) + t * sd(camp)/(sqrt(n))
  if(LL3 < mu & UL3 > mu) {
    sim3[i] = 1} else{
      sim3[i] = 0}
  # GM
  LL4 = gamma * exp(mean(log(camp)) - t * sd(log(camp))/(sqrt(n))) 
  UL4 = gamma * exp(mean(log(camp)) + t * sd(log(camp))/(sqrt(n))) 
  if(LL4 < mu & UL4 > mu){
    sim4[i] = 1} else{
      sim4[i] = 0}
}
out1 = mean(sim1)
out2 = mean(sim2)
out3 = mean(sim3)
out4 = mean(sim4)

# BETA(3,3)
a = 3
b = 3
mu =  a/(a+b) 
mug_x = function(x){
  out = exp(1/n * sum(log(x)))
  return(out)}
sim1 = matrix(NA, M)
sim2 = matrix(NA, M)
sim3 = matrix(NA, M)
sim4 = matrix(NA, M)
mux = matrix(NA, M)
mug = matrix(NA, M)
for(i in 1:M){
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c = n/(n-z)
  for(j in 1:M){
    camp = rbeta(n, a, b)  
    mux[j] = mean(camp)
    mug[j] = mug_x(camp)
  }
  mod1 = lm(mux ~ 0 + mug)
  gamma = coef(mod1)
  # Normal
  t = qt(1-alpha/2, n-1)
  LL1 = mean(camp) - t * sd(camp)/(sqrt(n))
  UL1 = mean(camp) + t * sd(camp)/(sqrt(n))
  if(LL1 < mu & UL1 > mu) {
    sim1[i] = 1 } else {
      sim1[i] = 0}
  # Bonett
  sigmab = sqrt(exp(log(c * var(camp))))
  LL2 = mean(camp) - z * sigmab/(sqrt(n))
  UL2 = mean(camp) + z * sigmab/(sqrt(n))
  if(LL2 < mu & UL2 > mu) {
    sim2[i] = 1} else{
      sim2[i] = 0}
  # Johnson
  mu3 = mean((camp - mean(camp))^3) 
  LL3 = (mean(camp) + mu3/(6*var(camp)*n)) - t * sd(camp)/(sqrt(n))
  UL3 = (mean(camp) + mu3/(6*var(camp)*n)) + t * sd(camp)/(sqrt(n))
  if(LL3 < mu & UL3 > mu) {
    sim3[i] = 1} else{
      sim3[i] = 0}
  # GM
  LL4 = gamma * exp(mean(log(camp)) - t * sd(log(camp))/(sqrt(n))) 
  UL4 = gamma * exp(mean(log(camp)) + t * sd(log(camp))/(sqrt(n))) 
  if(LL4 < mu & UL4 > mu){
    sim4[i] = 1} else{
      sim4[i] = 0}
}
out1 = mean(sim1)
out2 = mean(sim2)
out3 = mean(sim3)
out4 = mean(sim4)

###############################################################################
# TABELLA 7: INTERVALLI DI CONFIDENZA PER LA DIFFERENZA DI MEDIE
# LOGN(0,1)
n1 = 10
n2 = 15
M = 1000
mu = 0
mug_x = function(x, n){
  out = exp(1/n * sum(log(x)))
  return(out)}

sim1 = matrix(NA, M)
sim2 = matrix(NA, M)
sim3 = matrix(NA, M)
sim4 = matrix(NA, M)

mux1 = matrix(NA, M)
mux2 = matrix(NA, M)
mug1 = matrix(NA, M)
mug2 = matrix(NA, M)

for(i in 1:M){
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c1 = n1/(n1-z)
  c2 = n2/(n2-z)
  t = qt(1 - alpha/2, (n1+n2-2))
  for(j in 1:M){
    camp1 = rlnorm(n1, 0, 1)
    mux1[j] = mean(camp1)
    mug1[j] = mug_x(camp1, n1)}
  for(l in 1:M){
    camp2 = rlnorm(n2, 0, 1)
    mux2[l] = mean(camp2)
    mug2[l] = mug_x(camp2, n2)}
  
  mod1 = lm(mux1 ~ 0 + mug1)
  gamma1 = coef(mod1)
  mod2 = lm(mux2 ~ 0 + mug2)
  gamma2 = coef(mod2)
  
  # Normal
  LL1 = (mean(camp1)-mean(camp2)) - z * (sqrt((var(camp1)/n1) + (var(camp2)/n2)))
  UL1 = (mean(camp1)-mean(camp2)) + z * (sqrt((var(camp1)/n1) + (var(camp2)/n2)))
  if(LL1 < mu & UL1 > mu){
    sim1[i] = 1} else{
      sim1[i] = 0}
  amp1 = UL1-LL1
  # Bonett
  sigmab1 = (exp(log(c1 * var(camp1))))
  sigmab2 = (exp(log(c2 * var(camp2))))
  LL2 = (mean(camp1)-mean(camp2)) - z * (sqrt((sigmab1/n1) + (sigmab2/n2)))
  UL2 = (mean(camp1)-mean(camp2)) + z * (sqrt((sigmab1/n1) + (sigmab2/n2)))
  if(LL2 < mu & UL2 > mu){
    sim2[i] = 1} else{
      sim2[i] = 0}
  amp2=UL2-LL2
  # Jhonson
  mu13 = mean((camp1 - mean(camp1))^3)
  mu23 = mean((camp2 - mean(camp2))^3)
  a1 = mean(camp1) + (mu13/(6 * n1 * var(camp1)))
  a2 = mean(camp2) + (mu23/(6 * n2 * var(camp2)))
  LL3 = (a1 - a2) - t * sqrt((var(camp1)/n1 + var(camp2)/n2))
  UL3 = (a1 - a2) + t * sqrt((var(camp1)/n1 + var(camp2)/n2))
  if(LL3 < mu & UL3 > mu){
    sim3[i] = 1} else{
      sim3[i] = 0}
  amp3=UL3-LL3
  # GM
  LL = exp(mean(log(camp1)) - mean(log(camp2)) - t * (sqrt((var(log(camp1))/n1)) + (var(log(camp2))/n2)))
  UL = exp(mean(log(camp1)) - mean(log(camp2)) + t * (sqrt((var(log(camp1))/n1)) + (var(log(camp2))/n2)))
  LL4 = (gamma1/gamma2) * LL
  UL4 = (gamma1/gamma2) * UL
  if(LL4 < 1 & UL4 > 1){
    sim4[i] = 1} else{
      sim4[i] = 0}
  amp4=UL4-LL4
}
out1 = mean(sim1)
out2 = mean(sim2)
out3 = mean(sim3)
out4 = mean(sim4)


# UNIF(0,1)
n1 = 10
n2 = 15
M = 1000
mu = 0
mug_x = function(x, n){
  out = exp(1/n * sum(log(x)))
  return(out)}

sim1 = matrix(NA, M)
sim2 = matrix(NA, M)
sim3 = matrix(NA, M)
sim4 = matrix(NA, M)

mux1 = matrix(NA, M)
mux2 = matrix(NA, M)
mug1 = matrix(NA, M)
mug2 = matrix(NA, M)

for(i in 1:M){
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c1 = n1/(n1-z)
  c2 = n2/(n2-z)
  t = qt(1 - alpha/2, (n1+n2-2))
  for(j in 1:M){
    camp1 = runif(n1, 0, 1)
    mux1[j] = mean(camp1)
    mug1[j] = mug_x(camp1, n1)}
  for(l in 1:M){
    camp2 = runif(n2, 0, 1)
    mux2[l] = mean(camp2)
    mug2[l] = mug_x(camp2, n2)}
  
  mod1 = lm(mux1 ~ 0 + mug1)
  gamma1 = coef(mod1)
  mod2 = lm(mux2 ~ 0 + mug2)
  gamma2 = coef(mod2)
  
  # Normal
  LL1 = (mean(camp1)-mean(camp2)) - z * (sqrt((var(camp1)/n1) + (var(camp2)/n2)))
  UL1 = (mean(camp1)-mean(camp2)) + z * (sqrt((var(camp1)/n1) + (var(camp2)/n2)))
  if(LL1 < mu & UL1 > mu){
    sim1[i] = 1} else{
      sim1[i] = 0}
  amp1=UL1-LL1
  # Bonett
  sigmab1 = (exp(log(c1 * var(camp1))))
  sigmab2 = (exp(log(c2 * var(camp2))))
  LL2 = (mean(camp1)-mean(camp2)) - z * (sqrt((sigmab1/n1) + (sigmab2/n2)))
  UL2 = (mean(camp1)-mean(camp2)) + z * (sqrt((sigmab1/n1) + (sigmab2/n2)))
  if(LL2 < mu & UL2 > mu){
    sim2[i] = 1} else{
      sim2[i] = 0}
  amp2=UL2-LL2
  # Jhonson
  mu13 = mean((camp1 - mean(camp1))^3)
  mu23 = mean((camp2 - mean(camp2))^3)
  a1 = mean(camp1) + (mu13/(6 * n1 * var(camp1)))
  a2 = mean(camp2) + (mu23/(6 * n2 * var(camp2)))
  LL3 = (a1 - a2) - t * sqrt((var(camp1)/n1 + var(camp2)/n2))
  UL3 = (a1 - a2) + t * sqrt((var(camp1)/n1 + var(camp2)/n2))
  if(LL3 < mu & UL3 > mu){
    sim3[i] = 1} else{
      sim3[i] = 0}
  amp3 = UL3-LL3
  # GM
  LL = exp(mean(log(camp1)) - mean(log(camp2)) - t * (sqrt((var(log(camp1))/n1)) + (var(log(camp2))/n2)))
  UL = exp(mean(log(camp1)) - mean(log(camp2)) + t * (sqrt((var(log(camp1))/n1)) + (var(log(camp2))/n2)))
  LL4 = (gamma1/gamma2) * LL
  UL4 = (gamma1/gamma2) * UL
  if(LL4 < 1 & UL4 > 1){
    sim4[i] = 1} else{
      sim4[i] = 0}
  amp=UL4-LL4
}
out1 = mean(sim1)
out2 = mean(sim2)
out3 = mean(sim3)
out4 = mean(sim4)

# EXP(1)
n1 = 10
n2 = 15
M = 1000
mu = 0
mug_x = function(x, n){
  out = exp(1/n * sum(log(x)))
  return(out)}

sim1 = matrix(NA, M)
sim2 = matrix(NA, M)
sim3 = matrix(NA, M)
sim4 = matrix(NA, M)

mux1 = matrix(NA, M)
mux2 = matrix(NA, M)
mug1 = matrix(NA, M)
mug2 = matrix(NA, M)

for(i in 1:M){
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c1 = n1/(n1-z)
  c2 = n2/(n2-z)
  t = qt(1 - alpha/2, (n1+n2-2))
  for(j in 1:M){
    camp1 = rexp(n1, 1)
    mux1[j] = mean(camp1)
    mug1[j] = mug_x(camp1, n1)}
  for(l in 1:M){
    camp2 = rexp(n2, 1)
    mux2[l] = mean(camp2)
    mug2[l] = mug_x(camp2, n2)}
  
  mod1 = lm(mux1 ~ 0 + mug1)
  gamma1 = coef(mod1)
  mod2 = lm(mux2 ~ 0 + mug2)
  gamma2 = coef(mod2)
  
  # Normal
  LL1 = (mean(camp1)-mean(camp2)) - z * (sqrt((var(camp1)/n1) + (var(camp2)/n2)))
  UL1 = (mean(camp1)-mean(camp2)) + z * (sqrt((var(camp1)/n1) + (var(camp2)/n2)))
  if(LL1 < mu & UL1 > mu){
    sim1[i] = 1} else{
      sim1[i] = 0}
  amp1=UL1-LL1
  # Bonett
  sigmab1 = (exp(log(c1 * var(camp1))))
  sigmab2 = (exp(log(c2 * var(camp2))))
  LL2 = (mean(camp1)-mean(camp2)) - z * (sqrt((sigmab1/n1) + (sigmab2/n2)))
  UL2 = (mean(camp1)-mean(camp2)) + z * (sqrt((sigmab1/n1) + (sigmab2/n2)))
  if(LL2 < mu & UL2 > mu){
    sim2[i] = 1} else{
      sim2[i] = 0}
  amp2=UL2-LL2
  # Jhonson
  mu13 = mean((camp1 - mean(camp1))^3)
  mu23 = mean((camp2 - mean(camp2))^3)
  a1 = mean(camp1) + (mu13/(6 * n1 * var(camp1)))
  a2 = mean(camp2) + (mu23/(6 * n2 * var(camp2)))
  LL3 = (a1 - a2) - t * sqrt((var(camp1)/n1 + var(camp2)/n2))
  UL3 = (a1 - a2) + t * sqrt((var(camp1)/n1 + var(camp2)/n2))
  if(LL3 < mu & UL3 > mu){
    sim3[i] = 1} else{
      sim3[i] = 0}
  amp3=UL3-LL3
  # GM
  LL = exp(mean(log(camp1)) - mean(log(camp2)) - t * (sqrt((var(log(camp1))/n1)) + (var(log(camp2))/n2)))
  UL = exp(mean(log(camp1)) - mean(log(camp2)) + t * (sqrt((var(log(camp1))/n1)) + (var(log(camp2))/n2)))
  LL4 = (gamma1/gamma2) * LL
  UL4 = (gamma1/gamma2) * UL
  if(LL4 < 1 & UL4 > 1){
    sim4[i] = 1} else{
      sim4[i] = 0}
  amp4=UL4-LL4
}
out1 = mean(sim1)
out2 = mean(sim2)
out3 = mean(sim3)
out4 = mean(sim4)

# CHISQ(1)
n1 = 10
n2 = 15
M = 1000
mu = 0
mug_x = function(x, n){
  out = exp(1/n * sum(log(x)))
  return(out)}

sim1 = matrix(NA, M)
sim2 = matrix(NA, M)
sim3 = matrix(NA, M)
sim4 = matrix(NA, M)

mux1 = matrix(NA, M)
mux2 = matrix(NA, M)
mug1 = matrix(NA, M)
mug2 = matrix(NA, M)

for(i in 1:M){
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c1 = n1/(n1-z)
  c2 = n2/(n2-z)
  t = qt(1 - alpha/2, (n1+n2-2))
  for(j in 1:M){
    camp1 = rchisq(n1, 1)
    mux1[j] = mean(camp1)
    mug1[j] = mug_x(camp1, n1)}
  for(l in 1:M){
    camp2 = rchisq(n2, 1)
    mux2[l] = mean(camp2)
    mug2[l] = mug_x(camp2, n2)}
  
  mod1 = lm(mux1 ~ 0 + mug1)
  gamma1 = coef(mod1)
  mod2 = lm(mux2 ~ 0 + mug2)
  gamma2 = coef(mod2)
  
  # Normal
  LL1 = (mean(camp1)-mean(camp2)) - z * (sqrt((var(camp1)/n1) + (var(camp2)/n2)))
  UL1 = (mean(camp1)-mean(camp2)) + z * (sqrt((var(camp1)/n1) + (var(camp2)/n2)))
  if(LL1 < mu & UL1 > mu){
    sim1[i] = 1} else{
      sim1[i] = 0}
  amp1=UL1-LL1
  # Bonett
  sigmab1 = (exp(log(c1 * var(camp1))))
  sigmab2 = (exp(log(c2 * var(camp2))))
  LL2 = (mean(camp1)-mean(camp2)) - z * (sqrt((sigmab1/n1) + (sigmab2/n2)))
  UL2 = (mean(camp1)-mean(camp2)) + z * (sqrt((sigmab1/n1) + (sigmab2/n2)))
  if(LL2 < mu & UL2 > mu){
    sim2[i] = 1} else{
      sim2[i] = 0}
  amp2=UL2-LL2
  # Jhonson
  mu13 = mean((camp1 - mean(camp1))^3)
  mu23 = mean((camp2 - mean(camp2))^3)
  a1 = mean(camp1) + (mu13/(6 * n1 * var(camp1)))
  a2 = mean(camp2) + (mu23/(6 * n2 * var(camp2)))
  LL3 = (a1 - a2) - t * sqrt((var(camp1)/n1 + var(camp2)/n2))
  UL3 = (a1 - a2) + t * sqrt((var(camp1)/n1 + var(camp2)/n2))
  if(LL3 < mu & UL3 > mu){
    sim3[i] = 1} else{
      sim3[i] = 0}
  amp3=UL3-LL3
  # GM
  LL = exp(mean(log(camp1)) - mean(log(camp2)) - t * (sqrt((var(log(camp1))/n1)) + (var(log(camp2))/n2)))
  UL = exp(mean(log(camp1)) - mean(log(camp2)) + t * (sqrt((var(log(camp1))/n1)) + (var(log(camp2))/n2)))
  LL4 = (gamma1/gamma2) * LL
  UL4 = (gamma1/gamma2) * UL
  if(LL4 < 1 & UL4 > 1){
    sim4[i] = 1} else{
      sim4[i] = 0}
  amp4=UL4-LL4
}
out1 = mean(sim1)
out2 = mean(sim2)
out3 = mean(sim3)
out4 = mean(sim4)

# BETA(3,3)
n1 = 10
n2 = 15
a = 3
b = 3
mu =  a/(a+b) 
M = 1000
mu = 0
mug_x = function(x, n){
  out = exp(1/n * sum(log(x)))
  return(out)}

sim1 = matrix(NA, M)
sim2 = matrix(NA, M)
sim3 = matrix(NA, M)
sim4 = matrix(NA, M)

mux1 = matrix(NA, M)
mux2 = matrix(NA, M)
mug1 = matrix(NA, M)
mug2 = matrix(NA, M)

for(i in 1:M){
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c1 = n1/(n1-z)
  c2 = n2/(n2-z)
  t = qt(1 - alpha/2, (n1+n2-2))
  for(j in 1:M){
    camp1 = rbeta(n1, a, b)
    mux1[j] = mean(camp1)
    mug1[j] = mug_x(camp1, n1)}
  for(l in 1:M){
    camp2 = rbeta(n2, a, b)
    mux2[l] = mean(camp2)
    mug2[l] = mug_x(camp2, n2)}
  
  mod1 = lm(mux1 ~ 0 + mug1)
  gamma1 = coef(mod1)
  mod2 = lm(mux2 ~ 0 + mug2)
  gamma2 = coef(mod2)
  
  # Normal
  LL1 = (mean(camp1)-mean(camp2)) - z * (sqrt((var(camp1)/n1) + (var(camp2)/n2)))
  UL1 = (mean(camp1)-mean(camp2)) + z * (sqrt((var(camp1)/n1) + (var(camp2)/n2)))
  if(LL1 < mu & UL1 > mu){
    sim1[i] = 1} else{
      sim1[i] = 0}
  amp1=UL1-LL1
  # Bonett
  sigmab1 = (exp(log(c1 * var(camp1))))
  sigmab2 = (exp(log(c2 * var(camp2))))
  LL2 = (mean(camp1)-mean(camp2)) - z * (sqrt((sigmab1/n1) + (sigmab2/n2)))
  UL2 = (mean(camp1)-mean(camp2)) + z * (sqrt((sigmab1/n1) + (sigmab2/n2)))
  if(LL2 < mu & UL2 > mu){
    sim2[i] = 1} else{
      sim2[i] = 0}
  amp2=UL2-LL2
  # Jhonson
  mu13 = mean((camp1 - mean(camp1))^3)
  mu23 = mean((camp2 - mean(camp2))^3)
  a1 = mean(camp1) + (mu13/(6 * n1 * var(camp1)))
  a2 = mean(camp2) + (mu23/(6 * n2 * var(camp2)))
  LL3 = (a1 - a2) - t * sqrt((var(camp1)/n1 + var(camp2)/n2))
  UL3 = (a1 - a2) + t * sqrt((var(camp1)/n1 + var(camp2)/n2))
  if(LL3 < mu & UL3 > mu){
    sim3[i] = 1} else{
      sim3[i] = 0}
  amp3=UL3-LL3
  # GM
  LL = exp(mean(log(camp1)) - mean(log(camp2)) - t * (sqrt((var(log(camp1))/n1)) + (var(log(camp2))/n2)))
  UL = exp(mean(log(camp1)) - mean(log(camp2)) + t * (sqrt((var(log(camp1))/n1)) + (var(log(camp2))/n2)))
  LL4 = (gamma1/gamma2) * LL
  UL4 = (gamma1/gamma2) * UL
  if(LL4 < 1 & UL4 > 1){
    sim4[i] = 1} else{
      sim4[i] = 0}
  amp4=UL4-LL4
}
out1 = mean(sim1)
out2 = mean(sim2)
out3 = mean(sim3)
out4 = mean(sim4)

################################################################################
# GRAFICO 1
n = 50
M = 10000
g = n- 1
var_camp = numeric(M)
for (i in 1:M) {
  camp = rnorm(n, 0, 1)
  var_camp[i] = var(camp) * (n- 1) / n }
media_var = mean(var_camp)
varianza_var = var(var_camp)
hist(var_camp, breaks = 30, probability = TRUE,
     main = paste("Distribuzione Campionaria della Varianza"),
     xlab = "Varianza Campionaria", col = "lightblue",
     border = "black",
     ylab = "Densita")
# Curva teorica chi-quadro
curve(dchisq(x * (n- 1), g) * (n- 1),
      add = TRUE, col = "blue", lwd = 2)
text(1.4,1.5, "Curva teorica Chi-Quadro", col = "blue")

################################################################################
# GRAFICO 2
n = 50
M = 10000
g = n- 1
log_var = numeric(M)
for (i in 1:M) {
  camp = rnorm(n, 0, 1)
  var_camp = var(camp) * (n- 1) / n
  log_var[i] = log(var_camp)}
media_log = mean(log_var)
sd_log = sd(log_var)
hist(log_var, breaks = 30, probability = TRUE,
     main = paste("Distribuzione Campionaria del logaritmo
 della varianza"),
     xlab = "Log(Varianza Campionaria)", col = "lightblue",
     border = "black", ylab = "Densita")
# Curva normale teorica
curve(dnorm(x, media_log, sd_log),
      add = TRUE, col = "blue", lwd = 2)
text(0.4,1.5, "Curva teorica Normale", col = "blue")

################################################################################
# GRAFICI 3 E 4
n = 100
M = 10000
g = 1
var = 2 * g
smf = matrix(NA, M)
smn = matrix(NA, M)
lb0 = numeric(M)
ub0 = numeric(M)
lb1 = numeric(M)
ub1 = numeric(M)
for(i in 1:M){
  camp1 = rchisq(n, g)
  camp2 = rchisq(n, g)
  perc = 1/(2*(n-4)^(1/2))
  gamma41_1 = n*sum((camp1-mean(camp1))^4)/
    (sum((camp1-mean(camp1))^2))^2
  gamma41_2 = n*sum((camp2-mean(camp2))^4)/
    (sum((camp2-mean(camp2))^2))^2
  gamma42_1 = n*sum((camp1-mean(camp1,perc/2))^4)/
    (sum((camp1-mean(camp1))^2))^2
  gamma42_2 = n*sum((camp2-mean(camp2,perc/2))^4)/
    (sum((camp2-mean(camp2))^2))^2
  gamma43_1 = n*sum((camp1-median(camp1))^4)/
    (sum((camp1-mean(camp1))^2))^2
  gamma43_2 = n*sum((camp2-median(camp2))^4)/
    (sum((camp2-mean(camp2))^2))^2
  alpha = 0.05
  z = qnorm(1-alpha/2, 0, 1)
  c = n/(n-z)
  se11_1 = c * ((gamma41_1- (n-3)/n)/(n-1))^(1/2)
  se11_2 = c * ((gamma41_2- (n-3)/n)/(n-1))^(1/2)
  se12_1 = c * ((gamma42_1- (n-3)/n)/(n-1))^(1/2)
  se12_2 = c * ((gamma42_2- (n-3)/n)/(n-1))^(1/2)
  se13_1 = c * ((gamma43_1- (n-3)/n)/(n-1))^(1/2)
  se13_2 = c * ((gamma43_2- (n-3)/n)/(n-1))^(1/2)
  LLn = (var(camp2)/var(camp1)) * qf((alpha/2), n-1, n-1)
  lb0[i] = LLn
  ULn = (var(camp2)/var(camp1)) * qf((1-alpha/2), n-1, n-1)
  ub0[i] = ULn
  if(LLn < var & ULn > var) {smn[i] = 1} else{smn[i] = 0}
  f11 = sum((camp1- mean(camp1))^4)
  f12 = sum((camp2- mean(camp2))^4)
  mu4f1 = (f11 + f12)/(2*n)
  s2f1 = ((n-1) * var(camp1) + (n-1) * var(camp2))/(n+n)
  gamma4f1 = mu4f1/(s2f1^2)
  r1 = (2*n) / (gamma4f1- ((n-3)/(n-1)))
  r2 = (2*n) / (gamma4f1- ((n-3)/(n-1)))
  LLf1 = (var(camp2)/var(camp1)) * qf((alpha/2), r1, r2)
  lb1[i] = LLf1
  ULf1 = (var(camp2)/var(camp1)) * qf((1-alpha/2), r1, r2)
  ub1[i] = ULf1
  if(LLf1 < var & ULf1 > var) {smf[i] = 1} else {smf[i] = 0}}
outF1 = mean(smf)
outN = mean(smn)
# Grafico 3.1
data1 = data.frame(
  simu = 1:M,
  lower = lb0,
  upper = ub0,
  met = "NORMALE",
  cop = smn)
ggplot(data1, aes(x = simu, ymin = lower, ymax = upper)) +
  geom_linerange(aes(color = factor(cop))) +
  facet_wrap(~ met, ncol = 1) +
  geom_hline(yintercept = var, linetype = "dashed", color = "red") +
  labs(title = "Intervalli di confidenza NORMALE per la varianza",
       x = "Simulazioni", y = "Varianza") +
  scale_color_manual(values = c("1" = "blue", "0" = "lightblue"),
                     labels = c("Non Copre", "Copre")) +
  theme_minimal() +
  theme(legend.title = element_blank())
# Grafico 3.2
data2 = data.frame(
  simu = 1:M,
  lower = lb1,
  upper = ub1,
  met = "F",
  cop = smf)
ggplot(data2, aes(x = simu, ymin = lower, ymax = upper)) +
  geom_linerange(aes(color = factor(cop))) +
  facet_wrap(~ met, ncol = 1) +
  geom_hline(yintercept = var, linetype = "dashed", color = "red") +
  labs(title = "Intervalli di confidenza F per la varianza",
       x = "Simulazioni", y = "Varianza") +
  scale_color_manual(values = c("1" = "blue", "0" = "lightblue"),
                     labels = c("Non Copre", "Copre")) +
  theme_minimal() +
  theme(legend.title = element_blank())
