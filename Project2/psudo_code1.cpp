vec test_vec1 = zeros(col);
test_vec1(0) = 1;
vec test_vec2 = zeros(col);
test_vec2(1) = 1;

mat B = A;

while(True){
  transform(B)

  if ((times_looped%10 == 0) && test_enabled){
      if(dot(B*test_vec1,B*test_vec2) != 0){
          cout << "Orthogonality was not conserved!" << endl;
          exit(1);
      }
  }
}
