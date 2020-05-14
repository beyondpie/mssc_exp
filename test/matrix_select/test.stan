parameters{
    matrix[4, 3] s;
}

model{
    int myrows[5] = {1,1,1,1,1};
    int mycols[3] = {1,2,3};
    to_vector(s) ~ std_normal();
    print("sampling from std normal, ", s);
    /* print("matrix indexing: ", s[myrows, mycols]); */
    print("matrix indexing: ", s[:, myrows]);
}
