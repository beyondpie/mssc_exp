parameters{
    matrix[4, 3] s;
}

model{
    int myrows[5] = {1,1,1,1,1};
    int mycols[3] = {1,2,3};
    int mypos[2,2] = {{1,1}, {2,2}};
    to_vector(s) ~ std_normal();
    print("sampling from std normal, ", s);
    /* print("matrix indexing: ", s[myrows, mycols]); */
    print("matrix indexing: ", to_array_2d(s)[mypos]);
}
