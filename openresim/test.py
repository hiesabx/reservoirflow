import scipy.sparse as ss

# bug: value is not stored!
d = ss.lil_matrix((10, 1), dtype="double")
d[0, 0] += 1
print(d.toarray())  # [[0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]]

# bug: value is not stored!
d = ss.lil_matrix((1, 10), dtype="double")
d.data[0] = [1]
print(d.toarray())  # [[0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]]

# bug: value is not stored!
d = ss.lil_matrix((1, 10), dtype="double")
d[0] = 0
d.data[0] += [1]
print(d.toarray())

# bug: value is not stored!
d = ss.lil_matrix((1, 10), dtype="double")
d.data[0] += [2]
print(d.toarray())

# value is stored but not compatible
d = ss.lil_matrix((1, 10), dtype="double")
# d[0] = 1
d[0, 0] += 2
print(d.toarray())
