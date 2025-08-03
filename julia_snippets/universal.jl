using Oscar, QP

d = 5

F,zd = cyclotomic_field(d)

X = gpX(d)
Z = gpZ(d)
T = weil_U(ZN(d)[1 1; 0 1]) # not the real T-gate but modular T
S = T*Z^3 # Can pick S = T*X^a*Z^b for any (a,b) ≠ (0,0)
R = diagonal_matrix(F,vcat([1 for i in 1:d-1],[-1]))
H = weil_U(ZN(d)[0 1; -1 0]) # This is the modular S

C = matrix_group([S,H])
sl = matrix_group(H,T)
CR = matrix_group([S,H,R])
HR = matrix_group([H,R])
println((order(C), order(sl)))
#15000, 120
println((is_infinite(CR), is_infinite(HR)))

#quo(CR,HR) seems to run forever so maybe [CR:HR] = ∞ but should check their algorithm ?quo
