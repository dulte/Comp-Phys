max_val = N; //The larges element (known)

vec calc_max_index = largest_non_diag(M);
double calc_max = M(calc_max_index(0),calc_max_index(1));

if (calc_max != max_val){
    cout << "The function did not calculate the correct largest element!" << endl;
    exit(1);
}
