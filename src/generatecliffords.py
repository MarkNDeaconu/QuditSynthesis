from tests.qutrit_tests.qutrit import *

cliffords = set()

curr = H
ops = [H,S]
for i in range(100000):
    curr = random.choice(ops) * curr
    cliffords.add(curr)

print(len(cliffords))

with open('QuditSynthesis/cliffords3.pkl', 'wb') as f:
    pickle.dump(cliffords, f)

# with open('cliffords3.pkl', 'rb') as f:
#     cliffords = pickle.load(f)
#     cliffords = list(cliffords)


# final_cliffords = []

# for c in cliffords:
#     if not((-1) * c in final_cliffords):
#         c.string = ''
#         final_cliffords.append(c)


# with open('cliffords3.pkl', 'wb') as f:
#     pickle.dump(final_cliffords, f)

# print(len(final_cliffords))

