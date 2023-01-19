from math import comb

n = 34260 # number of T-cells we test
t = 1000000000 # total population of T-cell in drug product
m = 100000 # amount of T-cells that are RCL positive in drug product
p = m/t # probability of success (getting a T-cell with RCL in our RCL assay)
k = 3 # number of RCL T cells we need to get a positive result (at least 3 RCL T-cell which gives us 6 VSV-G copies)

prob = 0 # probability that we will get at least k RCL positive T-cells when there is m amout of RCL positive T-cells in p amount of total T-cells in drug product when we take out n amount of T-cells for the RCL assay

for i in range(k):
    prob += comb(n, i) * p**i * (1-p)**(n-i)

result = (1 - prob)
print(result)
