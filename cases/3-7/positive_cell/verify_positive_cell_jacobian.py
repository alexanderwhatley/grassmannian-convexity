import sympy as s
a=s.symbols('a0:12');word=[4,3,2,1,5,4,3,2,6,5,4,3]
T=s.eye(7)
for j,t in zip(word,a):
 for i in range(7):T[i,j-1]+=t*T[i,j]
Y=T[4:7,:];Z=Y[:,4:7];assert Z.det()==1
Q=s.simplify(Z.inv()*Y[:,:4]);flat=[Q[i,j] for i in range(3) for j in range(4)]
print('coordinates done',flush=True)
J=s.Matrix(flat).jacobian(a)
print('jacobian done',flush=True)
f=s.factor(J.det(method='domain-ge'))
print('Jdet=',f,flush=True)
assert f!=0
