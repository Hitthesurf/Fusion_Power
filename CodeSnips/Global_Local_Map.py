for j in range(d):
    F_Global[Map[j]] += F_Local[j]
    for i in range(d):
        A_Global_Sparse[Map[j], Map[i]] += A_Local[j, i]