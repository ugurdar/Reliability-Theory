rm = 0.2
rc = 0.3
rwc = 0.1
rbp = 0.4

# case 1
rf = rm * (1-(1-rwc*rbp)^2*(1-rwc))
r1 = 1- (1-rf)*(1-rc)
p1 = (1-rbp)*rbp
# case 2
r2 = r1
p2 =p1
# case 3
r3 = rm * (1-(1-rwc*rbp)^2)
p3 = (1-rbp)^2
# case 4
rff = rm *(1-(1-rwc*rbp)^2*(1-rwc)^2)
r4 = 1-(1-rff)*(1-rc)
p4 = rbp^2

rs = p1*r1+p2*r2+p3*r3+p4*r4
rs


a = rwc*rbp
first = 2*a - a^2
second = (2*rwc-rwc^2) * rc * (2*rbp * rbp^2)
t1  = rm * (first + second - first*second) 
t2 = rc * (2*rbp - rbp^2)  

t1 + t2 - t1 * t2


a = rwc*rbp  
first = a^4-2*a^3+3*a^2-2*a +1
second = a^4*rc^4 - 2*a^3*rc^3 + 3*a^2+rc^2 - 2*a*rc + 1
rm * (first*second + first*(1-second)+second*(1-first))
