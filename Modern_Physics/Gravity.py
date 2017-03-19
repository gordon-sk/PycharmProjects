d = .0745
r = .0465
s = .0173
M = 1.5
S = 2.83

CI = 8.8E-7
Ct = .072

w1 = .020966
w2 = .020942
lambda1 = .000956319
lambda2 = .00095694


trial_aprox = (d*pow(r,2)*s*pow(w2,2))/(4.0*M*S)

trial_better = (d*s*pow(r,2)*(pow(w2,2)+pow(lambda2,2))*(1+CI))/(4*M*S*(1-Ct))

print trial_better