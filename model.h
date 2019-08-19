/* Model */
/* y' = f(y,t) */

const int neqn = 2; //Number of equations

void f(double (&yt) [2], double y [], double t)
{
  yt[0] = 1.0;
  yt[1] = y[0] - y[1];
}


