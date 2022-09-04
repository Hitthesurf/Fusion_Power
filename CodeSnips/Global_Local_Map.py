for j in range(d):
    F_Global[Map[j]] += F_Local[j]
    for i in range(d):
        B_Global_Sparse[Map[j], Map[i]] += B_Local[j, i]