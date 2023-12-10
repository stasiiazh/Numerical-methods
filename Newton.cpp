double *Newton(double x1, double x2, int NIT, double E1, double E2) 
{
 int k = 1;
 cout << left << setw(4) << "k" << setw(20) << "d1" << setw(20) << "d2" << setw(20) << "x1" << setw(20) << "x2" << endl;
 double *F = new double[2];
 double **J = new double*[2];
 for (int i = 0; i < 2; i++)
 {
  J[i] = new double[2];
 }
 double *dX = new double[2];
 double x1k, x2k;
 double d1, d2;
 double tmp;
 do{
  F[0] = -f1(x1, x2);
  F[1] = -f2(x1, x2);
  Jcobi(x1, x2, J);
dX = Gauss(J, F, 2); 
  x1k = x1 + dX[0]; 
  x2k = x2 + dX[1]; 
  d1 = abs(f1(x1, x2)); 
  tmp = abs(f2(x1, x2)); 
  if (tmp > d1) 
  { 
   d1 = tmp; 
  } 
  d2 = abs(x1k - x1) / (x1k >= 1 ? x1k : 1); 
  tmp = abs(x2k - x2) / (x2k >= 1 ? x2k : 1); 
  if (tmp > d2) 
  { 
   d2 = tmp; 
  } 
  x1 = x1k; 
  x2 = x2k; 
  cout << left << setw(4) << k << setw(20) << d1 << setw(20) << d2 << setw(20) << x1 << setw(20) << x2 << endl; 
  if (k >= NIT) 
  { 
   cout << "IER=2\n"; 
   system("pause"); 
   exit(2); 
  } 
  k++; 
 } while (d1>E1 && d2>E2); 
 dX[0] = x1; 
 dX[1] = x2; 
 return dX ; 
} 
 
double *Newton_M(double x1, double x2, int NIT, double E1, double E2, double M) 
{ 
 int k = 1; 
 cout << left << setw(4) << "k" << setw(20) << "d1" << setw(20) << "d2" << setw(20) << "x1" << setw(20) << "x2" << endl; 
 double *F = new double[2]; 
 double **J = new double*[2]; 
 for (int i = 0; i < 2; i++) 
 { 
  J[i] = new double[2]; 
 } 
 double *dX = new double[2]; 
 double x1k, x2k; 
 double d1, d2; 
 double tmp; 
 do { 
  F[0] = -f1(x1, x2); 
  F[1] = -f2(x1, x2); 
  Jcobi_M(x1, x2, J, M); 
 
  dX = Gauss(J, F, 2); 
  x1k = x1 + dX[0]; 
  x2k = x2 + dX[1]; 
  d1 = max(abs(f1(x1, x2)), abs(f2(x1, x2))); 
  d2 = max(abs(x1k - x1) / max(1.0, abs(x1k)), abs(x2k - x2) / max(1.0, abs(x2k)));
  x1 = x1k; 
  x2 = x2k; 
  cout << left << setw(4) << k << setw(20) << d1 << setw(20) << d2 << setw(20) << x1 << setw(20) << x2 << endl; 
  if (k >= NIT) 
  { 
   cout << "IER=2\n"; 
   exit(2); 
  } 
  k++; 
 } while (d1 > E1 || d2 > E2); 
 dX[0] = x1; 
 dX[1] = x2; 
 return dX; 
}