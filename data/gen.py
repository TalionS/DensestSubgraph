import random
f = open('data/core2.txt','w')
n = 1000
m = 3000
f.write(str(n) + ' ' + str(m) + '\n')
D = {}
for t in range(m):
    a = random.randint(0,n - 1)
    b = random.randint(0,n - 1)
    while a == b or (a,b) in D or (b,a) in D:
        a = random.randint(0,n - 1)
        b = random.randint(0,n - 1)
    f.write(str(a) + ' ' + str(b) + '\n')
    D[(a,b)] = 1